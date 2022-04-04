def pipeline(func):
    def wrapper(*args, **kwargs):
        if args[0] is None:
            return None
        else:
            return func(*args, **kwargs)
    return wrapper


def userinterface(prompt: str):
    """
    A decorator used to make user-prompt-event-loops. Accepts a prompt for
    the user. Passes 'user_input' as the first positional argument to the
    wrapped function.

    E.g.,
    @userinterface
    >>> import cli
    >>> @cli.userinterface("give me some input: ")
    ... def foo(user_input, bar):
    ...     if user_input == bar:
    ...         return "That's a nice bar"
    ...     else:
    ...         raise ValueError("That's not as nice a bar")
    ...
    >>> foo("chocolate")
    --------------------------------------------------------------------------------
    give me some input: banana
    input error; raised the following:
    --------------------------------------------------------------------------------
    That's not as nice a bar
    give me some input: chocolate
    "That's a nice bar"
    >>> foo("never gonna give you up")
    --------------------------------------------------------------------------------
    give me some input: quit
    >>> None
    """
    quit_strings = {"quit", "exit"}

    def decorator(func):
        def wrapper(*args, **kwargs):
            print("-" * 80)
            while True:
                user_input = input(prompt)
                if user_input in quit_strings:
                    return None
                else:
                    try:
                        return func(user_input, *args, **kwargs)
                    except Exception as e:
                        print(f"input error; raised the following: ")
                        print(e)
                        print("-" * 80)

        return wrapper
    return decorator

