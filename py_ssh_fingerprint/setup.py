from distutils.core import setup, Extension

library_dirs = ['openbsd-compat/', '.', '/usr/lib']
include_dirs = ['.', 'openbsd-compat/', '/usr/include']
libraries = ['openbsd-compat', 'ssh', 'crypt', 'crypto']

module = Extension('sshkey',
        sources=['key.c', 'uuencode.c', 'buffer.c', 'xmalloc.c', 'log.c',
        'ssh-dss.c', 'compat.c', 'bufbn.c',  
        'sshkey.c'],
        library_dirs=library_dirs,
        include_dirs=include_dirs,
        libraries=libraries,
        extra_objects=['libssh.a', 'openbsd-compat/libopenbsd-compat.a'])

setup (name='sshkey',
        version='0.1',
        description='Python module for calling OpenSSH key_fingerprint()\
function',
        ext_modules=[module])
