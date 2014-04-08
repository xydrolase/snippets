/*
 * ssh.c
 *
 * Xin Yin <killkeeper@gmail.com>
 * 
 * A thin Python wrapper for the OpenSSH's key_fingerprint() function. 
 *
 */
#include <Python.h>

#include "key.h"

static PyObject *
sshkey_fingerprint(PyObject *self, PyObject *args)
{
    char *flat_key;
    char *fp_hex_digest;
    int fingerprint_type;
    Key *key;
    
    if (!PyArg_ParseTuple(args, "s|i", &flat_key, &fingerprint_type)){
        return NULL;
    }

    key = key_new(KEY_UNSPEC);
    if (key_read(key, &flat_key) != -1){
        fp_hex_digest = key_fingerprint(key, fingerprint_type, SSH_FP_HEX);
        key_free(key);
        
        return Py_BuildValue("s", fp_hex_digest);

    }
    else{
        return NULL;
    }
}

static PyMethodDef SshKeyMethods[] = {
    {"fingerprint", (PyCFunction)sshkey_fingerprint, METH_VARARGS,
        "Compute the fingerprint (hash hex digest) by given hash algorithm"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initsshkey(void){
    PyObject *module;
    module = Py_InitModule3("sshkey", SshKeyMethods,
            "Python module for calling OpenSSH key_fingerprint() function.");

    if (module == NULL){
        return;
    }

    PyModule_AddIntConstant(module, "SSH_FP_MD5", SSH_FP_MD5);
    PyModule_AddIntConstant(module, "SSH_FP_SHA1", SSH_FP_SHA1);
}
