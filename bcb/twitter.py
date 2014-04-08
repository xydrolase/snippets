#!/usr/bin/env python

import tweepy
from conf import OAUTH_CONSUMER_KEY, OAUTH_CONSUMER_SECRET

auth = tweepy.OAuthHandler(OAUTH_CONSUMER_KEY,
        OAUTH_CONSUMER_SECRET)

def publish_tweet(tweet):
    try:
        from conf import OAUTH_ACCESS_TOKEN_KEY, OAUTH_ACCESS_TOKEN_SECRET
    except ImportError:
        print 'Run python twitter.py to get your access token before publishing tweets!'
        return

    auth.set_access_token(
            OAUTH_ACCESS_TOKEN_KEY,
            OAUTH_ACCESS_TOKEN_SECRET
            )
    api = tweepy.API(auth)
    api.update_status(tweet)

def main():
    print """\
[ User authorization for BCB Reptile ]
--------------------------------------
Please access the following URL in your browser:
%s""" % auth.get_authorization_url()
    verifier = raw_input("Enter the verifier you got from Twitter:")
    if auth.get_access_token(verifier):
        fd = open("conf.py", "a")
        fd.write("""
## AUTOMATITALLY GENERATED BY BCB REPTILE
OAUTH_ACCESS_TOKEN_KEY = "%s"
OAUTH_ACCESS_TOKEN_SECRET = "%s"
""" % (auth.access_token.key,
    auth.access_token.secret))
        fd.close()

        print 'Successfully got your token for access.'

if __name__ == "__main__":
    main()