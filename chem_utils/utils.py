"""Helpful utility functions."""
from inspect import signature, Parameter
from typing import Any, Callable, Type

from docstring_parser import parse

from tap import Tap
from tap.utils import type_to_str


# TODO: move this to Tap and enforce argument order when dynamically creating a Tap class

def convert_to_tap(function: Callable) -> Type[Tap]:
    """Converts a function into a Tap class.

    :param function: The function whose parameters will be converted to a Tap class.
    :return: A Tap class with arguments from the function.
    """
    # Get signature from function
    sig = signature(function)

    # Parse function docstring
    docstring = parse(function.__doc__)

    # Get function description
    description = '\n'.join(filter(None, (docstring.short_description, docstring.long_description)))

    # Add arguments of function to the Tap class
    annotations, defaults = {}, {}

    for param_name, param in sig.parameters.items():
        # Get type of the argument
        if param.annotation != Parameter.empty:
            # Any type defaults to str (needed for dataclasses where all non-default attributes must have a type)
            if param.annotation is Any:
                annotations[param.name] = str
            # Otherwise, get the type of the argument
            else:
                annotations[param.name] = param.annotation

        # Get the default or required of the argument
        if param.default != Parameter.empty:
            defaults[param.name] = param.default

    # Set parameter docstrings in configure
    def configure(self):
        for param in docstring.params:
            req_status = f'default={defaults[param.arg_name]}' if param.arg_name in defaults else 'required'
            help = f'({type_to_str(annotations[param.arg_name])}, {req_status}) {param.description}'
            self.add_argument(f'--{param.arg_name}', help=help)

    # Create tap class
    tap_class = type(f'{function.__name__}_tap', (Tap,), {
        '__doc__': description,
        '__annotations__': annotations,
        'configure': configure,
    } | defaults)

    return tap_class
