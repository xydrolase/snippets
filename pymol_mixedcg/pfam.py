#!/usr/bin/env python
"""
    A python library for Pfam services

    Author: Xin Yin <xinyin at iastate dot edu>
"""

import urllib
import re
import threading
from time import sleep

class SearchThread(threading.Thread):
    class QueryError(Exception):
        error_codes = {
            502: 'Bad gateway. Problem occurred when scheduling or running\
the job',
            503: 'Service unavailable. Job accepted but on hold.',
            410: 'Gone. Queried job is deleted from search system.',
            500: 'Internal server error.'
        }
        def __init__(self, status):
            self.status = status
            self.value = error_codes.get(self.status, 'Unknown error.') 
        def __str__(self):
            return " : ".join([str(self.status), self.value])

    def __init__(sequence, callback, **kwargs):
        self.sequence = sequence
        self.callback = callback
        self.r_result_url = re.compile(
            r'<result_url>(.+)?</result_url>'
        )

        self.parameters = {'seq': sequence}
        for k, v in kwargs:
            self.parameters[k] = v

        threading.Thread.__init__(self)

    def run(self):
        f = urllib.urlopen(
            'http://pfam.sanger.ac.uk/search/sequence',
            urllib.urlencode(self.parameters)
        )

        if f.getcode() == 200:
            response = f.read()
            # Search for result_url, use regular expression in favor of XML
            # parsing, which is unnecessary.
            match = self.r_result_url.search(response)
            if match:
                result_url = match.group(1)

                while True:
                    fr = urllib.urlopen(result_url)
                    status_code = fr.getcode()
                    if status_code in [502, 503, 410, 500]:
                        raise QueryError(fr.getcode())
                    else status_code == 200:



                    sleep(2000) # interval for checking status


class Sequence:
    """Class that utilizes Pfam's sequence searching service"""

    class InvalidSequenceError(Exception):
        pass

    def __init__(self):
        self.rseq = re.compile(r'^[A-IK-NP-Z]+$', re.IGNORECASE)
        pass

    def search(self, sequence, callback):
        if isinstance(sequence, list):
            query_sequence = "".join(sequence).strip()
        elif isinstance(sequence, list):
            query_sequence = sequence.strip()

        if not self.rseq.match(query_sequence):
            raise InvalidSequenceError("Invalid protein sequence: %s" % 
                                       query_sequence)

        query = SearchThread(sequence, callback)
        query.start()

        
