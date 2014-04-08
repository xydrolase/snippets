#!/usr/bin/env python
""" BCB Reptile
     v0.1

     Xin Yin <xinyin@iastate.edu>

     A script that tracks updates of specified web pages,
     updating these changes to Twitter.
"""
import subprocess
import sys
import os
import time
import re
from twitter import publish_tweet
from urlfetch import fetch_courses

def strip_tags(html):
    """Remove all HTML tags"""
    plain_text = re.sub(r'</?\w[^>]*?>', '', html)
    return re.sub(r'&#\d+;', ' ', plain_text).strip()

class BcbReptile:
    def __init__(config):
        self.config = config
        self.repository_path = config['repository_path']
        self.cookies = {}

        self.re_diff_context = re.compile(r'^@@ .+ @@')

    def git_call(self, args): 
        """Wrapper for calling shell commands."""
        pipe = subprocess.Popen(
                args,
                shell=True,
                cwd=self.repository_path,
                stdout=subprocess.PIPE)
        
        _stdout = pipe.communicate()[0].strip()
        return _stdout, pipe.returncode

    def git_diffstat(self):
        """Look for diffs between current working tree (newly updated web pages)
        and the HEAD of index"""
        return self.git_call('git diff --numstat')

    def course_diff(self, course):
        """Find changes in webpage for given course.

        Only insertions are monitored. Modifications that took place in different
        locations will be treated as separate updates. 

        HTML tags will be stripped and the plain text updates will be returned.
        """
        diff, retcode = self.git_call(
                'git diff %s' % course
                )

        if retcode == 0 and diff:
            diff_lines = diff.split('\n')[4:]
            diff_context_indices = [lineno 
                    for lineno in xrange(len(diff_lines))
                    if self.re_diff_context.search(diff_lines[lineno])]

            if diff_context_indices[-1] != len(diff_lines):
                diff_context_indices.append(len(diff_lines))

            updated_contexts = []
            # Find locations in diff output where modifications take place.
            reduce(lambda x, y: updated_contexts.append(diff_body[x:y]) or y,
                    diff_context_indices) 

            for context in updated_contexts:
                inserted_lines = filter(lambda x: x.strip(),
                        [line[1:] for line in context
                            if line.startswith('+')]




    def follow_links(self, raw_html):
        pass



def follow_links(html):
    """Follow all HTML <a> tags that hint for document files. Synchronize the
    updated file to specified Dropbox repository"""
    match_iter = re.finditer(r'<a .*?href=([\'\"])(.+?)\1.*?>', html, 
            re.IGNORECASE)
    for match in match_iter:
        url = match.group(2)


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
            # Find locations in diff output where modifications take place.
            reduce(lambda x, y: chunks.append(diff_body[x:y]) or y,
                    chunks_pos) 
            
            # Separate diffs into several bins for multiple updates. 
            plain_text_diffs = []
            for diff_chunk in chunks:
                inserted_lines = filter(lambda x: x,
                        [line[1:] for line in diff_chunk 
                            if line.startswith('+')])
                
                #aggregated_changes = "".join(insert_lines)
                #follow_links(aggregated_changes)

                # Chop update message into 125 characters at most
                plain_changes = " ".join(map(strip_tags, inserted_lines))
                if len(plain_changes) > 125:
                    plain_changes = ''.join([plain_changes[:122], '...'])
                plain_text_diffs.append(plain_changes)

            return plain_text_diffs

def commit_all(path):
    """Once updates were tracked, we need to commit changes to webpages to the
    git repository for further monitor."""
    timestamp = time.strftime("%Y-%m-%d %H:%M")
    return git_call(
            'git commit -a -m "Update @ %s"' % timestamp,
            path)

def check_for_updates():
    """Look for any updates since our last synchornization."""
    path = os.path.join('/home/frog/bcb', 'watch')
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
    """Publish all updates to Twitter."""
    for course, changes in updates.iteritems():
        for update in changes:
            publish_tweet('[%s] %s' % (course, update))

def main():
    # Grab the latest web pages, and save them to working tree.
    fetch_courses()

    updates = check_for_updates()
    if updates:
        publish_updates(updates)

        git_path = os.path.join('/home/frog/bcb', 'watch')
        commit_all(git_path)

    log_path = os.path.join('/home/frog/bcb', 'output.log')
    log_file = open(log_path, 'a+')
    log = """[%s]
%s
-----
""" % (time.strftime("%Y-%m-%d %H:%M"), str(updates) if updates else 'No updates tracked.')
    log_file.write(log)
    log_file.close()

if __name__ == "__main__":
    main()
