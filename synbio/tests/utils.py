def raises(callable, arg_list, kwarg_list, exc_list):
    iterator = zip(arg_list, kwarg_list, exc_list)
    for (w, kw, Exc) in iterator:
        try:
            callable(*w, **kw)
            raise AssertionError(f"{callable.__name__} passed with args={w}, kwargs={kw}")

        except Exception as raised:
            if not isinstance(raised, Exc):
                raise AssertionError(f"{raised}")
