# coding=utf-8

"""
    stream.varint
    ~~~~~~~~~~~~~

    Encode and decode an integer up to 64 bit to/from 'Varint'. See Google
    Protobuf library documentation for more details about Varints.

    :copyright: (c) 2017 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

import sys

import click
from google.protobuf.internal.decoder import _DecodeVarint as decodeVarint
from google.protobuf.internal.encoder import _EncodeVarint as encodeVarint


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def cli():
    """Varint encoder/decoder."""
    pass


@cli.command('encode')
@click.argument('integer', nargs=1, type=int)
def cmd_encode(integer):
    """Encode an integer up to 64 bit to Varint."""
    encode(integer)


@cli.command('decode')
@click.argument('input_file', nargs=1, type=click.File('rb'))
def cmd_decode(input_file):
    """Decode an integer up to 64 bit from Varint."""
    decode(input_file)


def encode(value):
    """Output the encoded value to the standard output.

    Args:
        value (int):  the integer to be encoded.
    """
    encodeVarint(sys.stdout.buffer.write, value, True)


def decode(input_file):
    """Output the decoded value to the standard output.

    Args:
        input_file (file handler):  input file handler.
    """
    print(decodeVarint(input_file.read(), 0)[0])
