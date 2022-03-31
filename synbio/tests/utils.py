class TestNotImplemented(BaseException): pass

def raises(callable, args, kwargs):
    """
    A function that returns whatever exception is raised by a Callable. Returns None if no exception is raised

    Parameters
    ----------
        Callable callable: any callable (e.g., function, class, etc.)
        list args: list of positional arguments to pass to callable
        list kwargs: list of keyword arguments to pass to callable

    Returns
    -------
        Optional[Exception]: exception raised by callable (None if exits cleanly)
    """
    try:
        callable(*args, **kwargs)
        return None
    except Exception as raised:
        return raised
