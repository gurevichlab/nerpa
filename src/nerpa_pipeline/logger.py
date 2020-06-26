#!/usr/bin/env python

import sys

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += f'{s}\n'
        print(s)

    def warn(self, s):
        msg = f"WARNING: {s}\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = f"ERROR: {s}\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()

