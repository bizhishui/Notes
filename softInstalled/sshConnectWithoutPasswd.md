#Connect to remote server without enter passwd#
1. Generation keys 

ssh-keygen -t dsa

2. copy the key to remote server, only for the first time, if you want to connect to the server in another computer, just put the contents of local id_dsa.pub file to the last of authorized_keys file in the remote server

scp id_dsa.pub remote@address:$HOME/authorized_keys
