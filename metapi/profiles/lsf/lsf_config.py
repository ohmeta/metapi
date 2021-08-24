import shlex
from collections import OrderedDict
from itertools import chain
from typing import TextIO, Union, List, Any, Dict

import yaml


class Config:
    def __init__(self, data: Union[dict, None] = None):
        self._data = dict()
        if data is not None:
            for key, value in data.items():
                self._data[key] = self.concatenate_params(value)

    def __bool__(self) -> bool:
        return bool(self._data)

    def __contains__(self, item) -> bool:
        return item in self._data

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)

    @staticmethod
    def args_to_dict(args: str) -> Dict[str, str]:
        """Converts a string into a dictionary where key/value pairs are consecutive
        elements of the string.
        Eg '-J "2" -q 3' --> {'-J': '2', '-q': '3'}
        """
        args_iter = shlex.shlex(args, posix=True)
        args_iter.whitespace_split = True
        return OrderedDict(zip(args_iter, args_iter))

    @staticmethod
    def concatenate_params(params: Union[List[str], str]) -> str:
        if isinstance(params, str):
            return params
        return " ".join(filter(None, params))

    def default_params(self) -> str:
        return self.get("__default__", "")

    def params_for_rule(self, rulename: str) -> str:
        """Loads default + rule-specific arguments.
        Arguments specified for a rule override default-specified arguments.
        Shlex-joining is required to properly pass quoted escapes in yaml
        to the shell.
        """
        default_params = self.args_to_dict(self.default_params())
        rule_params = self.args_to_dict(self.get(rulename, ""))
        default_params.update(rule_params)
        return " ".join(map(shlex.quote, chain.from_iterable(default_params.items())))

    @staticmethod
    def from_stream(stream: TextIO) -> "Config":
        data = yaml.safe_load(stream)
        return Config(data)
