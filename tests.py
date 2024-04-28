import subprocess

key = bytes.fromhex('8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef')
input_data = bytes.fromhex('1122334455667700ffeeddccbbaa9988')

kuz = Kuznechik(key)

encrypted_data = kuz.encrypt(input_data)

decrypted_data = kuz.decrypt(encrypted_data)

if decrypted_data == input_data:
    print('Encryption and decryption successful')
else:
    print('Encryption and decryption failed')
