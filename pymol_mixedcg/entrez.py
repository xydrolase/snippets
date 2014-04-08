#!/usr/bin/env python

from Bio import Entrez

import sys
import urllib
import re
 
Entrez.email = "xinyin@iastate.edu"

PROTEIN_ENTRIES = ['1NES_E', '3ITI_A']
 
class EntrezRetriever:
    def __init__(self, db='gene'):
        self.db = db
        self.search_feature_header = re.compile(' {5}[A-za-z]+').search
        self.search_active_site_type = re.compile(r'/site_type="active').search
        self.search_active_sites = re.compile(r'order\((.+)\)').search
        self.cache_pool = {}

    def annotate(self, id):
        key = '/'.join((self.db, str(id)))
        if key in self.cache_pool:
            return self.cache_pool[key]
        else:
            raw = self.strip_tags(urllib.urlopen(
            'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val=%(id)s&db=\
%(db)s&retmode=html' % {'id' : id, 'db': self.db}).read()).splitlines(True)
            self.cache_pool[key] = raw
            return raw

    def search_term(self, term, max_results=10):
        handle = Entrez.esearch(self.db, retmax=max_results, term=term)
        try:
            record = Entrez.read(handle)
        except RuntimeError, e:
            print "An error occurred during searching database."
            print "ERROR: %s" % e
            sys.exit(1)

        if record['Count'] > 0:
            return record['IdList']
        else:
            return []

    def retrieve_sequence(self, id):
        raw_metadata = self.annotate(id)

        start_lineno = [lineno for lineno in xrange(len(raw_metadata))
                if raw_metadata[lineno].startswith("ORIGIN")][0]
        seqs = '\n'.join(filter(lambda x: x,
            map(lambda seq: ''.join(re.split(r'\s+', seq.strip())[1:]),
                raw_metadata[start_lineno+1:])))

        return seqs

    def retrieve_annotation(self, id_list):
        """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
        submit the data to NCBI) and esummary to retrieve the information. 
        Returns a list of dictionaries with the annotations."""
     
        request = Entrez.epost(self.db, id=",".join([str(id) for id in id_list]))
        try:
            result = Entrez.read(request)
        except RuntimeError, e:
            #FIXME: How generate NAs instead of causing an error with invalid IDs?
            print "An error occurred while retrieving the annotations."
            print "The error returned was %s" % e
            sys.exit(-1)
     
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        data = Entrez.esummary(db=self.db, webenv=webEnv, query_key = queryKey)
        annotations = Entrez.read(data)
     
        return annotations

    def retrieve_features(self, id, feature):
        raw_metadata = self.annotate(id)

        lines_with_headers = [lineno for lineno in xrange(len(raw_metadata))
                if self.search_feature_header(raw_metadata[lineno])]
        lines_with_features = filter(
            lambda lineno: raw_metadata[lineno].startswith('     %s' % feature),
            lines_with_headers)

        features_selected = []
        for site_lineno in lines_with_features:
            index = lines_with_headers.index(site_lineno)
            if index < len(raw_metadata) - 1:
                features_selected.append(raw_metadata[lines_with_headers[index]:\
                        lines_with_headers[index+1]])
            elif index == len(raw_metadata):
                features_selected(raw_metadata[lines_with_features[index]:])

        return features_selected

    def find_active_sites(self, all_sites, offset=0):
        # Only return the first one?
        for site in all_sites:
            if self.search_active_site_type(site[1]):
                return map(lambda acs:
                        int(acs)+offset,
                        self.search_active_sites(site[0]).group(1).split(','))

        return None

    def strip_tags(self, html):
        return re.sub('</?[^<>]+?>', '', html)

def main():
    if len(sys.argv) < 2:
        print 'Usage: %s [NCBI_PROTEIN_ENTRY] ...'
        sys.exit(1)

    entrezer = EntrezRetriever(db='protein')
    command = sys.argv[1]

    ruler = ' '.join(['  ^  '] * 3)
    for term in sys.argv[2:]:
        id_list = entrezer.search_term(term)
        if id_list:
            # Go after 1st term only
            print '>%s-%s' % (term, id_list[0])
            if command == 'sites':
                seq = ''.join(
                        entrezer.retrieve_sequence(id_list[0]).splitlines(1))
                active_sites = entrezer.find_active_sites(
                        entrezer.retrieve_features(id_list[0], 'Site'))

                print seq
                print active_sites
                print ','.join(seq[pos-3:pos+2] for pos in active_sites) 
                print ruler

            elif command in ('seq', 'sequence'):
                print entrezer.retrieve_sequence(id_list[0])

if __name__ == "__main__":
    main()
