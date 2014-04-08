#!/usr/bin/env python

import urllib2
import urllib
import os
import re

from conf import COURSES_TO_WATCH

TAGS_TO_WATCH = ('table', 'tr', 'td', 'th', # table
        'dd', 'dl', 'dd',
        'pre', 'code', 'backquote', 'p', 'div', 'h[1-6]', 
        'ul', 'li')

def fetch_courses():
    for course, uri in COURSES_TO_WATCH.iteritems():
        if type(uri) is str:
            response = urllib.urlopen(uri).read()
        elif type(uri) is dict:
            req = urllib2.Request(uri['uri'], headers=uri['headers'])
            f = urllib2.urlopen(req)
            if uri['follow'] and not f.geturl() == uri['uri']:
                # A redirect is required, now we need to parse the cookies...
                cookies_to_set = f.info().getheaders('Set-cookie')
                cookie_plain = ' '.join([ck[:ck.find(';')+1] 
                    for ck in cookies_to_set])
                
                # Initializing a new request
                req = urllib2.Request(
                        f.geturl(),
                        headers={'Cookie': cookie_plain}
                        )
                f = urllib2.urlopen(req)
                response = f.read()
            else:
                response = f.read()

        if response:
            pattern_end_tag = re.compile(r'</(%s)>(\S)' % 
                    ('|'.join(TAGS_TO_WATCH)), re.IGNORECASE)
            response = re.sub(pattern_end_tag,
                    r'</\1>\n\2', response)

            pattern_start_tag = re.compile(r'(\S)<(%s)' %
                    ('|'.join(TAGS_TO_WATCH)), re.IGNORECASE)
            response = re.sub(pattern_start_tag,
                    r'\1\n<\2', response)

            out_path = os.path.join('watch', "%s" % course)
            fd = open(out_path, 'w')
            fd.write(response)
            fd.close()
