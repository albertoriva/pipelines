#!/usr/bin/env python

# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import sys
import datetime

def timestamp():
    return datetime.datetime.now().isoformat()

class Logger():
  outfile = ""
  out = None
  echo = False

  def __init__(self, filename, overwrite=True, echo=False):
      self.setLogfile(filename, overwrite=overwrite)
      self.setEcho(echo)

  def close(self):
      if self.out:
          self.out.close()

  def setLogfile(self, filename, overwrite=True):
      self.outfile = filename
      mode = "w" if overwrite else "a"
      if self.out:
          self.out.close()
      if filename:
          self.out = open(filename, mode)
      else:
          self.out = None

  def setEcho(self, echo):
      if echo == 'stdout':
          self.echo = sys.stdout
      elif echo == 'stderr':
          self.echo = sys.stderr
      else:
          self.echo = None

  def log(self, message, *args):
      s = message.format(*args)
      if self.out:
          self.out.write("{}\t{}\n".format(timestamp(), s))
          self.out.flush()
          if self.echo:
              self.echo.write(s + "\n")
      return s

  def logStart(self, scriptname):
      self.log("#"*30)
      self.log("### Script {} started.", scriptname)

  def logEnd(self):
      self.log("### Script terminated.")
      self.log("#"*30)
