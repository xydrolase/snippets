#!/usr/bin/env python

import subprocess
import os
import sys
import re
from urlfetch import fetch_courses

def git_call(args, path):
    pipe = subprocess.Popen(
            args,
            shell=True,
            cwd=path,
            stdout=subprocess.PIPE)
    
    _stdout = pipe.communicate()[0].strip()
    return _stdout, pipe.returncode


def commit_all(path):
    git_call(
            'git init',
            path
            )
    git_call(
            'git add .',
            path
            )
    git_call(
            'git commit -a -m "Initial commit"',
            path
            )

def main():
    if not os.path.exists(os.path.join(os.getcwd(), 'watch')):
        os.mkdir('watch')

    fetch_courses()

    path = os.path.join(os.getcwd(), 'watch')
    commit_all(path)

if __name__ == "__main__":
    main()
