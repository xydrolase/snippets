#!/usr/bin/env python

from random import choice
from pwd import getpwnam
from subprocess import Popen, PIPE

import os

class UserExists(Exception):
    pass

def generate_password():
    alpha = ''.join([chr(x) for x in range(ord('A'), ord('Z') + 1)])
    digit = ''.join([str(x) for x in range(0, 10)])

    allowable_chars = '[]^#@' + alpha + alpha.lower() + digit
    return ''.join([choice(allowable_chars) for x in range(8)])

def system_invoke(args):
    return Popen(' '.join(map(str, args)), shell=True, cwd=os.getcwd(),
        stdout=PIPE).communicate()[0]

def useradd(user, password=None):
    try:
        user = getpwnam(user)
        raise UserExists
    except KeyError:
        if not password:
            password = generate_password()

        #system_invoke(["useradd", "--create-home", "--base-dir", "/home",
        # set bash as default shell
        system_invoke(["useradd", "--create-home", "-s", "/bin/bash", "--base-dir", "/home",
            "--password", password, user])

        return (user, password)
