from __future__ import division, print_function
import numpy as np


def check_equation(equation):
    if callable(equation):
        parameter = np.arange(5, dtype=np.float)
        answer = equation(parameter)
        if isinstance(answer, np.ndarray):
            if answer is None:
                return True
            if answer.size == parameter.size and answer.shape == parameter.shape:
                return True
    return False
