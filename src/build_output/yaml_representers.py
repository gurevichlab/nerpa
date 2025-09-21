from enum import Enum

def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)
