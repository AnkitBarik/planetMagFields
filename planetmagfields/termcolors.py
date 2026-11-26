#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal ANSI colouring for terminal messages.

Colour is used only when the target stream is an interactive terminal, so
output redirected to a file or pipe stays plain text. Set the environment
variable NO_COLOR (https://no-color.org) to disable colour, or FORCE_COLOR to
enable it even when the stream is not a terminal.
"""

import os
import sys

_CODES = {
    'red':    '31',
    'green':  '32',
    'yellow': '33',
    'cyan':   '36',
    'bold':   '1',
}


def use_color(stream=None):
    """True if ANSI colours should be written to stream (default stdout)."""
    if os.environ.get('NO_COLOR'):
        return False
    if os.environ.get('FORCE_COLOR'):
        return True
    stream = sys.stdout if stream is None else stream
    isatty = getattr(stream, 'isatty', None)
    try:
        return bool(isatty and isatty())
    except ValueError:      # closed stream
        return False


def colored(text, color, stream=None):
    """Return text wrapped in the ANSI code for color ('red', 'green',
    'yellow', 'cyan' or 'bold') if stream supports colour, else text."""
    if not use_color(stream):
        return str(text)
    return '\033[%sm%s\033[0m' % (_CODES[color], text)


def red(text, stream=None):
    return colored(text, 'red', stream)


def yellow(text, stream=None):
    return colored(text, 'yellow', stream)


def cyan(text, stream=None):
    return colored(text, 'cyan', stream)


def bold(text, stream=None):
    return colored(text, 'bold', stream)
