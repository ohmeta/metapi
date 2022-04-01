#!/usr/bin/env python3

import pandas as pd
import requests
import xmltodict
import argparse
from rich import print
from rich.console import Console

# https://github.com/Textualize/rich/issues/67
_console = Console()

class RichArgumentParser(argparse.ArgumentParser):
    def _print_message(self, message, file=None):
        _console.print(message)

    def add_argument_group(self, *args, **kwargs):
        group = super().add_argument_group(*args, **kwargs)
        group.title = f"[cyan]{group.title.title()}[/cyan]"
        return group


class RichRawTextHelpFormatter(argparse.RawTextHelpFormatter):
    def _split_lines(self, text, width):
        return [f"[yellow]{line}[/yellow]" for line in text.splitlines()]


# see: http://goo.gl/kTQMs
SYMBOLS = {
    'customary'     : ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'),
    'customary_ext' : ('byte', 'kilo', 'mega', 'giga', 'tera', 'peta', 'exa', 'zetta', 'iotta'),
    'iec'           : ('Bi', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi', 'Yi'),
    'iec_ext'       : ('byte', 'kibi', 'mebi', 'gibi', 'tebi', 'pebi', 'exbi', 'zebi', 'yobi'),
}


def bytes2human(n, format='%(value).1f %(symbol)s', symbols='customary'):
    n = int(n)
    if n < 0:
        raise ValueError("n < 0")
    symbols = SYMBOLS[symbols]
    prefix = {}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i+1)*10
    for symbol in reversed(symbols[1:]):
        if n >= prefix[symbol]:
            value = float(n) / prefix[symbol]
            return format % locals()
    return format % dict(symbol=symbols[0], value=n)


def human2bytes(s):
    init = s
    num = ""
    while s and s[0:1].isdigit() or s[0:1] == '.':
        num += s[0]
        s = s[1:]
    if num != "":
        num = float(num)
    else:
        raise ValueError(f"can't covert {s} to float")
    letter = s.strip()
    #print(letter)
    for name, sset in SYMBOLS.items():
        if letter in sset:
            break
    else:
        if (letter == 'k') or (letter == "m") or (letter == "g"):
            # treat 'k' as an alias for 'K' as per: http://goo.gl/kTQMs
            sset = SYMBOLS['customary']
            letter = letter.upper()
        else:
            raise ValueError("can't interpret %r" % init)
    prefix = {sset[0]:1}
    for i, s in enumerate(sset[1:]):
        prefix[s] = 1 << (i+1)*10
    return int(num * prefix[letter])


def generate_xml(http_link):
    print(f'''Parsing: {http_link}\n''')
    r = requests.get(http_link)
    if "xml" in r.headers['content-type']:
        print(f'''Success: get XML document from the link: {http_link}\n''')
        return r.text
    else:
        print(f'''Error: can't get XML document from the link: {http_link}\nExiting\n''')
        return None


def estimate_size(xml_str, output=None):
    xml_dict = xmltodict.parse(xml_str)
    if "ListBucketResult" in xml_dict:
        file_info_df = pd.DataFrame(xml_dict["ListBucketResult"]["Contents"])\
            .astype({"Size": int})\
            .sort_values(["Size", "Key"])
        print(file_info_df)

        if output is not None:
            file_info_df.to_csv(output, sep="\t", index=False)

        total_size = sum(file_info_df["Size"])
        total_size_human = bytes2human(total_size)
        print(f'''\nTotal file size is: {total_size}''')
        print(f'''\nTotal file size is: {total_size_human}''')
    else:
        print("\nError: parse XML document error\nExiting")


def main():
    parser = RichArgumentParser("Estimate T2T data size")
    parser.add_argument("--http-link", dest="http_link",
                        default="https://s3-us-west-2.amazonaws.com/human-pangenomics?/delimiter=/&prefix=T2T",
                        help="T2T file/directory S3 link, default:\nhttps://s3-us-west-2.amazonaws.com/human-pangenomics?/delimiter=/&prefix=T2T")
    parser.add_argument("--output", dest="output", default=None,
                        help="a tsv file contains file information, default=None")
    args = parser.parse_args()

    xml_str = generate_xml(args.http_link)
    estimate_size(xml_str, args.output)


if __name__ == "__main__":
    main()
