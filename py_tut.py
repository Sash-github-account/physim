#!/usr/bin/env python
import sys
import optparse
import random
import time
from user import User
import urllib
import subprocess

vldpwd = 'mypass'

ippwd = raw_input('enter psswd: ')
counter = 0

parser = optparse.OptionParser()
parser.add_option('-n', '--name', dest='name', help='your name')
parser.add_option('-a', '--age', dest='age', help='your age', type=int)

(options, args) = parser.parse_args()

def game_func():
    ch_case=raw_input('enter what to do : 1. num game 2.str list manip')

    if ch_case == '1':
        answer = random.randint(1,10)
        num = 0

        while num != answer:
            num = int(raw_input('what number am i thinking of? '))

            if(num != answer):
                print "wrong"

        print "correct"

    elif ch_case == '2':
        struser = raw_input('enter a list of user names separated by commas')
        arr=struser.split(',')

        for user in arr:
            trimuser = user.strip()
            trimuserr = user.rstrip()
            trimuserl = user.lstrip()

            firstinit = trimuser[:1]
            lastinit = trimuser[1:2]
            lastname = trimuser[1:]

            print 'user: ' + user
            print 'ltrim: '+ trimuserl
            print 'rtrim: '+ trimuserr
            print 'first init: '+ firstinit
            print 'last init: '+ lastinit
            print 'last name: '+ lastname

    elif ch_case == '3':
        user_op_choice = '2'
        theuser = []
        while user_op_choice != 'q':
            user_op_choice = raw_input('1.display 2. load an entry 3. save 4. load info q.quit')
            if user_op_choice != 'q':
                if(user_op_choice == '2'):
                    u = User()
                    u.loadfromIP()
                    theuser.append(u)

                elif (user_op_choice == '1'):
                    for u in theuser:
                        u.display()

                elif (user_op_choice == '3'):
                    for u in theuser:
                        u.save()

                elif (user_op_choice == '4'):
                    for u in theuser:
                        u.loadfromfile()

    elif ch_case == '4':
        properties = {}

        properties['protocol'] = 'http'
        properties['host'] = 'www.google.com'
        properties['port'] = '80'
        properties['path'] = '/trends/'

        url = properties['protocol'] + '://' + \
              properties['host'] + ':' + \
              properties['port'] + \
              properties['path']

        print 'reading url', url
        response = urllib.urlopen(url)
        print response.read()

    elif ch_case == '5':
        subprocess.call("ls -l", shell = True)

        proc = subprocess.Popen(['tail', '-10', 'STDcellLib_map.py'], stdout=subprocess.PIPE)

        for line in proc.stdout.readlines():
            print line.rstrip()


while 1:
    time.sleep(1)
    counter = counter + 1
    if ippwd == vldpwd :
        print 'access'
        if options.name is None:
            options.name = raw_input('enter name : ')

        if options.age is None:
            options.age = int(raw_input('enter age : '))

        sayhello= 'hello' + options.name +','

        if options.age == 100:
            sayage = 'u r 100 yrs'
        elif options.age < 100:
            sayage = 'ull b 100 in '+ str(100 - options.age) + 'yrs'
        else:
            sayage = 'u turned 100 '+ str(options.age - 100)+ ' yrs ago'

        game = raw_input('need somethng?(y/n):')
        if(game == 'y'):
            game_func()
        else:
            print 'script has been running for', counter, 'seconds...'
            sys.exit(0)

    else :
        print 'access denied'
        print 'script has been running for', counter, 'seconds...'
        sys.exit(0)
print sayhello, sayage
print 'script has been running for', counter, 'seconds...'
sys.exit(0)
