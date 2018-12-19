def has_equal_elements(left, right, epsilon=None):
    """
    Returns true whenever two iterable objects represent the same set of elements in the same order
    Actual type of the lists does not matter
    """
    __tracebackhide__ = True
    if len(right) != len(left):
        return False
    if epsilon is not None:
        return all([abs(left_elem - right_elem) <= epsilon for left_elem, right_elem in zip(left, right)])
    return all([left_elem == right_elem for left_elem, right_elem in zip(left, right)])
