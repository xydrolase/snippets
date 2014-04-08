#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
deploy_wordpress.py

Quickly deploy a current version WordPress blog system onto current server 
for designated user.

Author: xydrolase

"""

from argparse import ArgumentParser
from config.account import mysql

from unix import useradd, system_invoke, UserExists

import sys
import os

import MySQLdb

nginx_template = """
server {{
    listen 80;
    server_name {server_name};
    index index.html index.htm index.php;
    root /var/www/{user_name}/;

    client_max_body_size 4M;
    client_body_buffer_size 128k;

    if (-d $request_filename) {{
        rewrite ^/(.*)([^/])$ http://$host/$1$2/ permanent;
    }}

    location ~ .*\.php?$ {{
        fastcgi_pass unix:/var/run/php5-fpm.sock;
        fastcgi_index index.php;
        include fastcgi_params;
    }}
}}
"""

wordpress_template = """
<?php
define('DB_NAME', '{db_name}');
define('DB_USER', '{user_name}');
define('DB_PASSWORD', '{password}');
define('DB_HOST', 'localhost');
define('DB_CHARSET', 'utf8');
define('DB_COLLATE', '');

{secret_key}

$table_prefix  = 'wp_';

define('WPLANG', '');
define('WP_DEBUG', false);

if ( !defined('ABSPATH') )
	define('ABSPATH', dirname(__FILE__) . '/');

require_once(ABSPATH . 'wp-settings.php');
"""

def info(category, msg):
    print >>sys.stderr, "[{0}] {1}".format(category, msg)

def setup_latest_wordpress(path, db_name, user_name, password,
        rename="wordpress"):
    # unpack tarball
    info("WORDPRESS", "Downloding latest version of WordPress.")

    import urllib2
    response = urllib2.urlopen("http://wordpress.org/latest.tar.gz")
    raw_stream = response.read()

    if not os.path.exists(os.path.join(path, rename)):
        tarball_path = os.path.join(path, 'wp_latest.tar.gz')
        tarball = open(tarball_path, 'wb+')
        tarball.write(raw_stream)
        tarball.close()

        os.chdir(path)

        info("WORDPRESS", "Decompressing tar ball.")

        system_invoke(["tar", "zxf", tarball_path])
        if not rename == "wordpress":
            system_invoke(["mv", "wordpress", rename])

        os.remove("wp_latest.tar.gz")

    if not os.path.exists(os.path.join(path, rename, 'wp-config.php')):
        # generate wordpress secret-key
        info("WORDPRESS", "Generating secret key.")

        salt = urllib2.urlopen("https://api.wordpress.org/secret-key/1.1/salt/")
        salt_config = salt.read()

        # modify wordpress config
        info("WORDPRESS", "Initializing WordPress configurations.")

        os.chdir(os.path.join(path, rename))
        wp_conf = wordpress_template.format(db_name=db_name, 
                user_name=user_name,
                password=password,
                secret_key=salt_config)
        wp_file = open("wp-config.php", "w+")
        wp_file.write(wp_conf)
        wp_file.close()
    else:
        info("WORDPRESS", "WordPress already installed.")

def create_www(user_name):
    info("WWW", "Creating web directory for {0}".format(user_name))
    try:
        os.mkdir("/var/www/{0}".format(user_name))
    except OSError:
        return

    system_invoke(["chown", "-R", "{0}:{0}".format(user_name), 
        "/var/www/{0}".format(user_name)])

def init_nginx_conf(server_name, user_name):
    nginx_conf = nginx_template.format(server_name=server_name,
            user_name=user_name)
    nginx_dir = "/home/{0}/sites-available".format(user_name)

    try:
        info('NGINX', "Creating directory {0}".format(nginx_dir))
        os.mkdir(nginx_dir, 0770)
    except OSError:
        if os.path.exists(os.path.join(nginx_dir, user_name)):
            return

    conf_path = os.path.join(nginx_dir, user_name)
    info("NGINX", "Setting up user specific nginx configs {0}".format(
        conf_path))
    conf_file = open(conf_path, 'w+')
    conf_file.write(nginx_conf)
    conf_file.close()

    # grant directory ownership
    info("NGINX", "Enabling user specific nginx configs.")

    system_invoke(["chown", "-R", "{0}:{0}".format(user_name), nginx_dir])
    # create symbolic link for user specific nginx config
    if not os.path.exists("/etc/nginx/sites-enabled/{0}".format(
        user_name)):
        system_invoke(["ln", "-s", conf_path, \
                "/etc/nginx/sites-enabled/{0}".format(user_name)])
    # reload nginx configurations
    system_invoke(["service", "nginx", "reload"])

    return True

def main():
    parser = ArgumentParser(
            description="Automatic Wordpress deployment script")
    parser.add_argument("-b", dest="blog_name", default="wordpress",
            help="Directory name of which WordPress is installed.")
    parser.add_argument("-u", "--user", dest="user_name", required=True,
            help="Specify the user to which Wordpress is deployed.")
    parser.add_argument("-p", "--password", dest="password", default=None,
            help="Set user password (generate random password if omitted).")
    parser.add_argument("-s", "--server-name", dest="server_name", 
            default=None, 
            help="The domain to bind to designated user" + \
                    " (use <username>.monoer.me by default).")

    args = parser.parse_args()

    if args.server_name is None:
        args.server_name = "{0}.monoer.me".format(args.user_name)

    # (1) create system login
    try:
        info("USERADD", "Adding login to system: {0}".format(args.user_name))
        (login, passwd) = useradd(args.user_name, args.password)
        info("USERADD", "Login {0} created with password   {1}".format(
            login, passwd))
        args.password = passwd

    except UserExists:
        info("USERADD", "Login {0} already exists".format(args.user_name))

    # (2) creating database
    db = MySQLdb.connect(host=mysql['host'], 
            user=mysql['username'],
            passwd=mysql['password'])

    cursor = db.cursor()

    cursor.execute("FLUSH PRIVILEGES;")
    try:
        cursor.execute("CREATE USER '{0}'@'{1}' IDENTIFIED BY '{2}';".format(
            args.user_name, mysql['host'], args.password))
    except:
        pass

    # grant usage permission
    sql_grant = "GRANT USAGE ON *.* TO '{0}'@'{1}' IDENTIFIED BY '{2}'" + \
            "WITH MAX_QUERIES_PER_HOUR 0 MAX_CONNECTIONS_PER_HOUR 0 " + \
            "MAX_UPDATES_PER_HOUR 0 MAX_USER_CONNECTIONS 0;"

    cursor.execute(sql_grant.format(args.user_name, 
        mysql['host'], args.password))

    # gran all permission on the user's own table, add by loading
    sql_grant = "GRANT ALL PRIVILEGES ON `{0}` . * TO '{0}'@'localhost' " + \
            "WITH GRANT OPTION;"
    cursor.execute(sql_grant.format(args.user_name))

    # create a database with the same name as the user
    cursor.execute("CREATE DATABASE IF NOT EXISTS `{0}`;".format(
        args.user_name))

    cursor.close()
    db.close()

    # (3) installing wordpress
    create_www(args.user_name)
    setup_latest_wordpress("/var/www/{0}".format(args.user_name),
        args.user_name, args.user_name, args.password, 
        rename=args.blog_name)

    # (4) initialize nginx configuration
    init_nginx_conf(args.server_name, args.user_name)

if __name__ == "__main__":
    main()
