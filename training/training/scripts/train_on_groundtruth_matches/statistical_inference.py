from typing import Tuple
import numpy as np

LinFunc = Tuple[float, float]  # Define LinFunc type alias

NUM_POINTS = 1000  # Number of points for discretization
UNIFORM = np.ones(NUM_POINTS) / NUM_POINTS  # Uniform distribution


def validate_linear_function(f: LinFunc):
    """
    Validates that the linear function f(X) = a * X + b maps [0, 1] to [0, 1].

    Args:
        f: Linear function f(X) = a * X + b as a tuple (a, b).

    Raises:
        ValueError: If the function does not map [0, 1] to [0, 1].
    """
    a, b = f
    f_0 = a * 0 + b  # f(0)
    f_1 = a * 1 + b  # f(1)

    if not (0 <= f_0 <= 1) or not (0 <= f_1 <= 1):
        raise ValueError(
            f"The function f(X) = {a} * X + {b} does not map [0, 1] to [0, 1]. "
            f"f(0) = {f_0}, f(1) = {f_1}."
        )


def transform_distribution(prior: np.ndarray, f: LinFunc) -> np.ndarray:
    """
    Transforms a prior distribution using the linear function f(X) = a * X + b.

    Args:
        prior: The prior distribution for X.
        f: Linear function f(X) = a * X + b as a tuple (a, b).

    Returns:
        Transformed prior distribution for f(X).
    """
    # Validate the linear function
    #validate_linear_function(f)

    a, b = f
    num_points = len(prior)
    transformed_prior = np.zeros_like(prior)

    for i, prob in enumerate(prior):
        x = i / num_points  # Implied grid point for X
        fx = a * x + b
        fx_index = int(fx * num_points)
        transformed_prior[fx_index] += prob / abs(a)

    # Normalize the transformed prior
    transformed_prior /= transformed_prior.sum()
    return transformed_prior


def inverse_function(f: LinFunc) -> LinFunc:
    """
    Returns the inverse of a linear function f(X) = a * X + b.

    Args:
        f: Linear function f(X) = a * X + b as a tuple (a, b).

    Returns:
        Inverse function f^-1(Y) = (Y - b) / a as a tuple (a_inv, b_inv).
    """
    a, b = f
    return 1 / a, -b / a


def expectation(prior: np.ndarray) -> float:
    """
    Computes the expectation of a discrete distribution.

    Args:
        prior: The discrete distribution.

    Returns:
        Expectation of the distribution.
    """
    return float(np.sum(prior * np.linspace(0, 1, len(prior))))


def perform_bayesian_update(prior: np.ndarray, failures: int, successes: int) -> np.ndarray:
    """
    Performs Bayesian inference using Bernoulli likelihood.

    Args:
        prior: The prior distribution.
        successes: Number of successes in the data.
        failures: Number of failures in the data.

    Returns:
        Posterior distribution after Bayesian update.
    """
    num_points = len(prior)
    likelihood = np.linspace(0, 1, num_points) ** successes * (1 - np.linspace(0, 1, num_points)) ** failures
    unnormalized_posterior = prior * likelihood
    return unnormalized_posterior / unnormalized_posterior.sum()


def bayesian_inference(
        prior: np.ndarray,
        f: LinFunc,
        failures: int,
        successes: int
) -> np.ndarray:
    """
    Generalized Bayesian inference for an RV X with a linear transformation f(X) = a * X + b.

    Args:
        prior: The prior probability distribution for X.
        f: Linear function f(X) = a * X + b as a tuple (a, b).
        successes: Number of successes in the Bernoulli data.
        failures: Number of failures in the Bernoulli data.

    Returns:
        Posterior distribution for X.
    """
    # Step 1: Transform prior for X into prior for f(X)
    prior_fx = transform_distribution(prior, f)

    # Step 2: Perform Bayesian inference for f(X)
    posterior_fx = perform_bayesian_update(prior_fx, failures, successes)

    # Step 3: Transform posterior for f(X) back to posterior for X
    f_inv = inverse_function(f)
    posterior = transform_distribution(posterior_fx, f_inv)

    return posterior


if __name__ == "__main__":
    # Example Usage:
    num_points = 1000
    prior = np.ones(num_points) / num_points  # Uniform prior

    # Linear transformation for P(event) = 1 - A * (1 - p)
    A = 0.898
    f = (A, 1 - A)  # f(p) = A * p + (1 - A)

    # Successes and failures
    successes = 3
    failures = 1

    # Perform Bayesian inference
    posterior = bayesian_inference_general(prior, f, successes, failures)

    # Compute posterior mean
    posterior_mean = np.sum(posterior * np.linspace(0, 1, num_points))
    print(f"Posterior mean for PKS presence probability: {posterior_mean:.4f}")
