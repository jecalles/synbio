from synbio.codes import Code


def codesavvy(func):
    def wrapper(*args, **kwargs):
        # handle code input
        code = kwargs.get("code", None)
        if code is None:
            kwargs["code"] = Code()
        elif isinstance(code, dict):
            kwargs["code"] = Code(code)
        else:
            raise TypeError("code must be a dict or dict-like obj")
        return func(*args, **kwargs)
    return wrapper
