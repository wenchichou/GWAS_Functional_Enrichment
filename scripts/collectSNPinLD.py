#!/usr/bin/python
import sys
import urllib
import urllib2
from bs4 import BeautifulSoup

# Test the html link to get the result of the LD block table
# http://archive.broadinstitute.org/mammals/haploreg/haploreg.php?query=rs733381&ldThresh=0.8&ldPop=EUR
url = 'http://archive.broadinstitute.org/mammals/haploreg/haploreg.php'
values = { 'query': sys.argv[1],
		   'ldThresh' : '0.8',
		   'ldPop' : 'EUR'}

data = urllib.urlencode(values)
req = urllib2.Request(url, data)
response = urllib2.urlopen(req)
the_page = response.read()
soup = BeautifulSoup(the_page, 'html.parser')

outfile = open(sys.argv[2], 'w')
for row in soup.findAll('table')[3].findAll('tr')[1:]:
	variantLink = row.findAll('td')[4].text
	outfile.write("%s\n" % variantLink)

outfile.close()

