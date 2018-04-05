# this should be run as python 2
import sys
if sys.version_info[0] > 2:
    sys.exit("Must be using Python 2 -- Exitting")
    
#!/usr/bin/python
import string
import hashlib
from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file
def get_pw(directory, file):
    """
    Collect password string for BRENDA from directory/file
    """
    with open(directory+file, 'r') as f:
        for line in f:
            l = line.rstrip()
            return l    

        # input parameters
database_directory = '/home/atarzia/psp/sequence_db/test_dataset/'

import io

output = io.StringIO()
output.write(u'http://www.brenda-enzymes.org/soap/brenda.wsdl')

# Retrieve file contents -- this will be
# u'First line.\nSecond line.\n'
contents = output.getvalue()
print contents

# ID

username = 'andrew.tarzia@adelaide.edu.au'
pw_dir = '/home/atarzia/psp/brenda_details/'
pw_file = 'pw.txt'

# URLS
wsdl = contents
url = "http://www.brenda-enzymes.org/soap/brenda_server.php"
        

#1) Usage with WSDL (for extracting the URL of the endpoint)
password = get_pw(pw_dir, pw_file)
client = WSDL.Proxy(wsdl)
parameters = username+","+password+",ecNumber*1.1.1.1#organism*Homo sapiens#"
resultString = client.getKmValue(parameters)
print resultString

