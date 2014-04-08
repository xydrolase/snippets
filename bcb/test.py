#!/usr/bin/env python

import subprocess
import sys
import os
import urllib
from conf import *
import time
import re
from twitter import publish_tweet

def git_call(args, path):
    pipe = subprocess.Popen(
            args,
            shell=True,
            cwd=path,
            stdout=subprocess.PIPE)
    
    _stdout = pipe.communicate()[0].strip()
    return _stdout, pipe.returncode

def diff(path):
    return git_call(
            'git diff --numstat',
            path)

def strip_tags(html):
    plain_text = re.sub(r'<[^>]*?>', '', html)
    return re.sub(r'&#\d+;', ' ', plain_text).strip()

def diff_course(path, file):
    diff, retcode =  git_call(
            'git diff ' + file,
            path)

    if retcode == 0:
        if diff:
            diff_body = diff.split('\n')[4:]
            re_diff_chunk = re.compile(r'^@@ .+ @@')
            chunks_pos = [idx for idx in range(len(diff_body))
                    if re_diff_chunk.search(diff_body[idx])]

            if chunks_pos[-1] != len(diff_body):
                chunks_pos.append(len(diff_body))
            
            chunks = []
            reduce(lambda x, y: chunks.append(diff_body[x:y]) or y,
                    chunks_pos) 
            
            plain_text_diffs = []
            for diff_chunk in chunks:
                inserted_lines = filter(lambda x: x,
                        [line[1:] for line in diff_chunk 
                            if line.startswith('+')])
                
                plain_changes = strip_tags(" ".join(inserted_lines))[:120]
                plain_text_diffs.append(plain_changes)

            return plain_text_diffs

def commit_all(path):
    timestamp = time.strftime("%Y-%m-%d %H:%M")
    return git_call(
            'git commit -a -m "Update @ %s"' % timestamp,
            path)

def check_for_updates():
    path = os.path.join(os.getcwd(), 'watch')
    diffs, retcode = diff(path)
    if not retcode == 0:
        return None
    
    if diffs:
        courses_updated = [c.split('\t')[-1]
                for c in diffs.split('\n')]
        
        course_diffs = {}
        for course in courses_updated:
            changes = diff_course(path, course)
            if changes:
                course_diffs[course] = changes

        return course_diffs

def publish_updates(updates):
    for course, changes in updates.iteritems():
        for update in changes:
            publish_tweet('[%s] %s' % (course, update))

def main():
    updates = check_for_updates()
    if updates:
        print updates
        publish_updates(updates)

if __name__ == "__main__":
    main()
