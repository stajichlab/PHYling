#!/usr/bin/env python3
import sys

from html.parser import HTMLParser
import re
import urllib.request

class JGIHTMLParser(HTMLParser):
    """JGI HTML Parser"""
    myreadytable = 0
    myreadyaref = 0
    genomepref = []
    genomename = {}
    mytags = []
    jgiurl = re.compile('https://genome\.jgi\.doe\.gov/(\S+)')

#    def __init__(self):
#        self.myready = 0
#        self.genomes = []
#        self.mytages = []
#        self.jgiurl = re.compile('https://genomes\.jgi\.doe\.gov/(\S+)')


    def handle_starttag(self, tag, attrs):
        if tag == "table":
            for attr in attrs:
                if attr[0] == "class":
                    self.myreadytable = 1
        elif self.myreadytable == 1 and tag == "a":
            for attr in attrs:
                m = self.jgiurl.search(attr[1])
                if m:
                    self.genomepref.append(m.group(1))
                    self.myreadyaref = 1
                else:
                    self.myreadyaref = 0
                    
        self.mytags.append( tag )

    def handle_endtag(self, tag):
        if tag == "table":
            self.myreadytable = 0
        elif tag == "a":
            self.myreadyaref = 0
        JGIHTMLParser.mytags.pop()

    def handle_data(self, data):
        if self.myreadyaref == 1:
            self.genomename[self.genomepref[-1]] = data
            dataup = re.sub(" ","_",data)
            print(self.genomepref[-1], dataup)


if len(sys.argv) > 1:
    htmlfile = sys.argv[1]
else:
    print("Need an input file")
    exit()

parser = JGIHTMLParser()

with open(htmlfile,"rt") as fh:
    data = []
    for line in fh:
        data.append(line.strip())
    datastr = "".join(data)
    parser.feed(datastr)
