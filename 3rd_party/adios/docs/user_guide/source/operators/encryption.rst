**********
Encryption
**********

The Encryption Operator uses the :ref:`Plugins` interface.
This operator uses `libsodium <https://doc.libsodium.org/>`_ for encrypting and decrypting data.
If ADIOS can find libsodium at configure time, this plugin will be built.

This operator will generate a secret key and encrypts the data with the key and a nonce as described in the libsodium `secret key cryptography docs <https://doc.libsodium.org/secret-key_cryptography/secretbox>`_.
The key is saved to the specified ``SecretKeyFile`` and will be used for decryption. The key should be kept confidential since it is used to both encrypt and decrypt the data.

Parameters to use with the Encryption operator:

============================== ===================== ===========================================================
 **Key**                       **Value Format**       **Explanation**
============================== ===================== ===========================================================
 PluginName                     string                Required. Name to refer to plugin, e.g., ``MyOperator``
 PluginLibrary                  string                Required. Name of shared library, ``EncryptionOperator``
 SecretKeyFile                  string                Required. Path to secret key file
============================== ===================== ===========================================================
