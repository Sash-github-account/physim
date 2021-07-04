#!/usr/bin/env python
import sys

class User:
      name = ""
      age = 0
      height = 0
      weight = 0

      def display(self):
      	  print 'User info'
	  print 'user name: ' + self.name
	  print 'user age: '+ str(self.age)
	  print 'user height: ' + str(self.height)
	  print 'user weight: ' + str(self.weight)

      def loadfromIP(self):
            self.name = raw_input('enter name: ')
            self.age = int(raw_input('enter age: '))
            self.height = int(raw_input('enter height: '))
            self.weight = int(raw_input('enter weight: '))

      def save(self):
            f = open('user.info', 'w')
            f.write(self.name + '\n')
            f.write(str(self.age)+ '\n')
            f.write(str(self.height) + '\n')
            f.write(str(self.weight) + '\n')
            f.close()

      def loadfromfile(self):
            f = open('user.info', 'r')
            self.name = f.readline().rstrip()
            self.age = int(f.readline())
            self.height = float(f.readline())
            self.weight = float(f.readline())
