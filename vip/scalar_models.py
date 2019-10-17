import abc
import numpy as np

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

tf.enable_v2_behavior()
tfd = tfp.distributions


class ScalarModel(abc.ABC):
    """An abstract base class for Scalar Models.

    Samples are laid out as particles x variables, which we will call
    z-shape.

    * p is the target probability.
    * q is the distribution that we're fitting to p.

    * grad_z is laid out as (particle, variable, parameter)
    """

    def __init__(self, initial_params, variable_count, particle_count):
        assert initial_params.ndim == 1
        self.q_params = np.full((variable_count, len(initial_params)), initial_params)
        self.particle_count = particle_count
        # The current stored sample.
        self.z = None
        # The gradient of z with respect to the parameters of q.
        self.grad_z = None
        # The stochastic gradient of log sum q for z.
        self.grad_sum_log_q = None

    @property
    def variable_count(self):
        return self.q_params.shape[0]

    @property
    def param_count(self):
        return self.q_params.shape[1]

    def clear_sample(self):
        self.z = None
        self.grad_z = None
        self.grad_sum_log_q = None

    def suggested_step_size(self):
        return np.average(np.abs(self.q_params), axis=0) / 100

    def elbo_estimate(self, log_p, particle_count=None):
        """A naive Monte Carlo estimate of the ELBO.

        log_p must take an argument of z-shape.
        """
        if particle_count is None:
            particle_count = self.particle_count
        z, log_prob = self.sample(particle_count, log_prob=True)
        return (np.sum(log_p(z) - log_prob)) / particle_count

    @staticmethod
    def _chain_rule(grad_log_p_z, grad_z):
        return (
            np.tensordot(grad_log_p_z.transpose(), grad_z, axes=1)
            .diagonal()
            .transpose()
        )

    @staticmethod
    def _slow_chain_rule(grad_log_p_z, grad_z):
        particle_count, variable_count, param_count = grad_z.shape
        result = np.zeros((variable_count, param_count))
        for variable in range(variable_count):
            for param in range(param_count):
                for particle in range(particle_count):
                    result[variable, param] += (
                        grad_log_p_z[particle, variable]
                        * grad_z[particle, variable, param]
                    )
        return result

    def elbo_gradient_using_current_sample(self, grad_log_p_z, test_chain_rule=False):
        assert self.grad_z is not None
        if not np.all(np.isfinite(grad_log_p_z)):
            raise Exception(
                "Infinite gradient given to elbo_gradient_using_current_sample."
            )
        if test_chain_rule:
            if not np.allclose(
                self._chain_rule(grad_log_p_z, self.grad_z),
                self._slow_chain_rule(grad_log_p_z, self.grad_z),
            ):
                print("Warning: chain rule isn't close. Here's the difference:")
                print(
                    self._chain_rule(grad_log_p_z, self.grad_z)
                    - self._slow_chain_rule(grad_log_p_z, self.grad_z)
                )
        unnormalized_result = (
            self._chain_rule(grad_log_p_z, self.grad_z) - self.grad_sum_log_q
        )
        return unnormalized_result / self.particle_count

    @abc.abstractmethod
    def mode_match(self, modes):
        pass

    @abc.abstractmethod
    def sample(self, particle_count, log_prob=False):
        pass

    @abc.abstractmethod
    def sample_and_prep_gradients(self):
        pass


class LogNormalModel(ScalarModel):
    """A log-normal model with hand-computed gradients."""

    def __init__(self, initial_params, variable_count, particle_count):
        super().__init__(initial_params, variable_count, particle_count)
        self.name = "LogNormal"

    @property
    def mu(self):
        return self.q_params[:, 0]

    @property
    def sigma(self):
        return self.q_params[:, 1]

    def mode_match(self, modes):
        """Some crazy heuristics for mode matching with the given branch
        lengths."""
        log_modes = np.log(np.clip(modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(modes, 1e-6, 1 - 1e-6))
        # TODO do we want to have the lognormal variance be in log space? Or do we want
        # to clip it?
        self.q_params[:, 1] = -0.1 * biclipped_log_modes
        self.q_params[:, 0] = np.square(self.sigma) + log_modes

    def sample(self, particle_count, log_prob=False):
        z = np.random.lognormal(
            self.mu, self.sigma, (particle_count, self.variable_count)
        )
        if log_prob:
            return z, self.log_prob(z)
        else:
            return z

    def log_prob(self, z):
        """Return a log probability for each of the particles given in z."""
        log_z = np.log(z)
        # Here's the fancy stable version.
        # ratio = np.exp(np.log((np.log(np.log(z)-mu)**2) - np.log(sigma**2))
        ratio = (log_z - self.mu) ** 2 / (2 * self.sigma ** 2)
        # Below we're summing over axis 1, which is the variable axis.
        result = (
            -0.5 * np.log(2 * np.pi)
            - np.sum(log_z, axis=1)
            - np.sum(np.log(self.sigma))
            - np.sum(ratio, axis=1)
        )
        return result

    def sample_and_prep_gradients(self):
        self.z = self.sample(self.particle_count)
        epsilon = (np.log(self.z) - self.mu) / self.sigma
        self.grad_z = np.empty((self.particle_count, self.variable_count, 2))
        self.grad_z[:, :, 0] = self.z
        self.grad_z[:, :, 1] = self.z * epsilon
        self.grad_sum_log_q = np.empty((self.variable_count, 2))
        # To get the gradients below, recall that log z is mu+sigma*epsilon in the
        # reparametrization trick. If we substitute that into the log density of the
        # normal distribution and take the derivatives we get this.
        self.grad_sum_log_q[:, 0] = -1.0 * self.particle_count
        self.grad_sum_log_q[:, 1] = (
            -np.sum(epsilon, axis=0) - self.particle_count / self.sigma
        )


def exponential_factory(params):
    return tfd.Exponential(rate=params[:, 0])


def gamma_factory(params):
    return tfd.Gamma(concentration=tf.exp(params[:, 0]), rate=tf.exp(params[:, 1]))


def lognormal_factory(params):
    return tfd.LogNormal(loc=params[:, 0], scale=params[:, 1])


def truncated_lognormal_factory(params):
    exp_shift = tfp.bijectors.Chain(
        [
            tfp.bijectors.AffineScalar(shift=-tf.math.exp(params[:, 2])),
            tfp.bijectors.Exp(),
        ]
    )
    return tfd.TransformedDistribution(
        distribution=tfd.TruncatedNormal(
            loc=params[:, 0], scale=params[:, 1], low=params[:, 2], high=999
        ),
        bijector=exp_shift,
        name="TruncatedLogNormal",
    )


class TFScalarModel(ScalarModel):
    """An object to model a collection of scalar variables with a given
    distribution type via TensorFlow.

    See ScalarModel for more information.
    """

    def __init__(self, q_factory, initial_params, variable_count, particle_count):
        super().__init__(initial_params, variable_count, particle_count)
        self.q_factory = q_factory
        self.name = "TF" + self.q_factory(np.array([initial_params])).name

    def mode_match(self, modes):
        """Some crazy heuristics for mode matching with the given branch
        lengths."""
        log_modes = np.log(np.clip(modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(modes, 1e-6, 1 - 1e-6))
        # TODO do we want to have the lognormal variance be in log space? Or do we want
        # to clip it?
        if self.name == "TFLogNormal":
            self.q_params[:, 1] = -0.1 * biclipped_log_modes
            self.q_params[:, 0] = np.square(self.q_params[:, 1]) + log_modes
        elif self.name == "TFTruncatedLogNormal":
            self.q_params[:, 1] = -0.1 * biclipped_log_modes
            self.q_params[:, 0] = np.square(self.q_params[:, 1]) + log_modes
            self.q_params[:, 2] = -5
        elif self.name == "TFGamma":
            self.q_params[:, 1] = np.log(-60.0 * biclipped_log_modes)
            self.q_params[:, 0] = np.log(1 + modes * self.q_params[:, 1])
        else:
            print("Mode matching not implemented for " + self.name)

    def sample(self, particle_count, log_prob=False):
        q_distribution = self.q_factory(self.q_params)
        z = q_distribution.sample(particle_count).numpy()
        if log_prob:
            return z, np.sum(q_distribution.log_prob(z), axis=1)
        else:
            return z

    def sample_and_prep_gradients(self):
        """Take a sample from q and prepare a gradient of the sample and of log
        q with respect to the parameters of q."""
        with tf.GradientTape(persistent=True) as g:
            tf_params = tf.constant(self.q_params, dtype=tf.float32)
            g.watch(tf_params)
            q_distribution = self.q_factory(tf_params)
            tf_z = q_distribution.sample(self.particle_count)
            q_term = tf.math.reduce_sum(q_distribution.log_prob(tf_z))
        self.z = tf_z.numpy()
        # The Jacobian is laid out as particles x edges x edges x params.
        self.grad_z = np.sum(g.jacobian(tf_z, tf_params).numpy(), axis=2)
        self.grad_sum_log_q = g.gradient(q_term, tf_params).numpy()
        del g  # Should happen anyway but being explicit to remember.
        return self.z


def of_name(name, *, variable_count, particle_count):
    def build_tf_model(q_factory, initial_params):
        return TFScalarModel(
            q_factory, np.array(initial_params), variable_count, particle_count
        )

    if name == "lognormal":
        return LogNormalModel(np.array([-2.0, 0.5]), variable_count, particle_count)
    if name == "tf_lognormal":
        return build_tf_model(lognormal_factory, [-2.0, 0.5])
    if name == "tf_truncated_lognormal":
        return build_tf_model(truncated_lognormal_factory, [-1.0, 0.5, 0.1])
    if name == "tf_gamma":
        return build_tf_model(gamma_factory, [1.3, 3.0])

    raise Exception(f"Model {name} not known.")