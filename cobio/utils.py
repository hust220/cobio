import math
import requests, collections
import os, sys
import tarfile
from threading import Lock
from concurrent.futures import ThreadPoolExecutor
import time
import subprocess
from io import BytesIO, StringIO
import networkx as nx
import matplotlib.pyplot as plt
import json
import networkx.algorithms.isomorphism as iso

lock = Lock()

def distance(p1, p2):
  dx = p1[0] - p2[0]
  dy = p1[1] - p2[1]
  dz = p1[2] - p2[2]
  return math.sqrt(dx * dx + dy * dy + dz * dz)

def tar_add(tar_file, files, verbose=False):
  tar = tarfile.open(tar_file, mode='a')
#  f = open(new_file)

  try:
    filenames = ' '.join(files)
    if verbose:
      print(f'Adding {filenames} to {tar_file} ... ')
    for new_file in files:
      tar.add(new_file)
#    tarinfo = tarfile.TarInfo(name=new_file)
#    tarinfo.size = len(f.read())
#    tarinfo.mtime = time.time()
#    tar.addfile(tarinfo, fileobj=f)
  except Exception as e:
    print(f'Failed!')
    print(str(e))

#  f.close()
  tar.close()

def tar_list(filename):
  t = tarfile.open(filename, mode='r')
  ls = []
  for m in t.getmembers():
    ls.append(m.name)
  t.close()
  return ls

def file_to_list(filename):
  f = open(filename)
  lines = f.read().splitlines()
  f.close()
  return lines

class Compressor:
  def __init__(self, verbose=False, max_files=100):
    self.to_compress = 0
    self.max_files = max_files
    self.tasks = dict()
    self.todo = dict()
    self.verbose = verbose

  def compress(self, tasks):
    for tar in tasks:
      filenames = ' '.join(tasks[tar])
      if self.verbose:
        print(f'[*] Compiling {filenames} to {tar} ...')
      tar_add(tar, tasks[tar], verbose=self.verbose)
      for filename in tasks[tar]:
        os.remove(filename)

  def submit(self, tar_file, files):
    with lock:
      if isinstance(files, str):
#        print('Is string')
        files = [files]
      self.to_compress += len(files)
#      print(self.to_compres)
      if self.verbose:
        print(f'[**] {self.to_compress} files to compress ...')
      if tar_file not in self.tasks:
        self.tasks[tar_file] = []
      self.tasks[tar_file].extend(files)
      if self.to_compress > self.max_files:
        self.compress(self.tasks)
        self.to_compress = 0
        self.tasks = dict()

  def __del__(self):
    if self.to_compress > 0:
      self.compress(self.tasks)

def run_tasks(task, *wrong_args, ilist='', open_files=False, itar='', otar='', skip=lambda name,files:False, ext='', args=dict(), verbose=False, nthreads=1):
  """
  Run Task for each file in a Tar file or each name in a list file.

  Parameters
  task: The task function
  ilist: The input list file that contains a string at each line
  itar: The input tar file.
  otar: The output tar file.
  """

  if len(wrong_args) > 0:
    raise Exception('Error: Wrong arguments of run_tasks function!')

  # get existing_ss
  if os.path.exists(otar):
    if verbose:
      print(f'Checking {otar} ...')
    existing_files = set(tar_list(otar))
  else:
    existing_files = set()

  if itar != '':
    if not os.path.exists(itar):
      raise Exception(f'Error: {itar} not exists!')

    if verbose:
      print(f'[*] Opening {itar} ...')
    t = tarfile.open(itar, mode='r')
    args['compressor'] = Compressor(verbose=verbose)
    with ThreadPoolExecutor(nthreads) as executor:
      if verbose:
        print(f'[*] Getting members of {itar} ...')
      members = t.getmembers()
      itask = 1
      for m in members:
  #      name = os.path.splitext(os.path.basename(m.name))[0]
        if ext != '':
          name = os.path.splitext(os.path.basename(m.name))[0]
          if f'{name}{ext}' in existing_files:
            continue
        elif skip(m.name, existing_files):
          continue
        f = t.extractfile(m)
        content = f.read()
        if verbose:
          print(f'[*] Submit task for {m.name}')
        executor.submit(task, m.name, content, args)
        itask += 1
    t.close()
    if verbose:
      print(f'[*] All done.')

  elif ilist != '':
    if not os.path.exists(ilist):
      raise Exception(f'Error: {ilist} not exists!')

    f = open(ilist, mode='r')
    args['compressor'] = Compressor(verbose=verbose)
    with ThreadPoolExecutor(nthreads) as executor:
      itask = 1
      for line in f:
        line = line.strip()
        if ext != '':
          name = os.path.splitext(os.path.basename(line))[0]
          if f'{name}{ext}' in existing_files:
            continue
        elif skip(line, existing_files):
          continue
        content = open(line).read() if open_files else ''
        if verbose:
          print(f'[*] Submit task for {line}')
        executor.submit(task, line, content, args)
        itask += 1
    f.close()
    if verbose:
      print(f'[*] All done.')

class R2Converter:
  pair_symbols = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]

  def determine_level(self, ss, a, b):
      if ss[a] != '.' or ss[b] != '.':
          return -1
      level = 0
      while True:
          score = 0
          for i in range(a + 1, b):
              if ss[i] == self.pair_symbols[level][0]:
                  score += 1
              elif ss[i] == self.pair_symbols[level][1]:
                  score -= 1
          if score == 0:
              return level
          level += 1
          if level >= len(self.pair_symbols):
              return -1

  def bp2db(self, pairs, n):
    ss = ['.' for i in range(n)]

    for i,j in pairs:
      if i > j:
        i,j = j,i

      level = self.determine_level(ss, i, j)

      if level != -1:
        ss[i] = self.pair_symbols[level][0]
        ss[j] = self.pair_symbols[level][1]

    return ''.join(ss)



def mc2db(io):
  converter = R2Converter()
  residues_mark = 'Residue conformations'
  basepairs_mark = 'Base-pairs'

  residues = dict()
  basepairs = dict()
  iresidue = 0
  pairs = []
  context = 'other'
  for line in io:
    line = line.strip()
    if line.endswith('--------'):
      if line.startswith(residues_mark):
        context = 'residues'
      elif line.startswith(basepairs_mark):
        context = 'basepairs'
      else:
        context = 'other'
    else:
      if context == 'residues':
        residue = line.strip().split()[0]
        residues[residue] = iresidue
        iresidue += 1
      elif context == 'basepairs':
        if 'Ww/Ww' in line:
          r1, r2 = line.strip().split()[0].split('-')
          i1 = residues[r1]
          i2 = residues[r2]
          pairs.append((i1,i2))
  return converter.bp2db(pairs, iresidue)

class R2Graph:
  pair_symbols = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]
  chain_symbols = ['&', ' ']

  def __init__(self):
    self.open_symbols = [p[0] for p in self.pair_symbols]
    self.close_symbols = [p[1] for p in self.pair_symbols]


  def get_open_symbol(self, close_symbol):
    ind = self.close_symbols.index(close_symbol)
    return self.open_symbols[ind]

  def find_helices(self, db):
    helices = []
    n = len(db)
    open_stacks_dict = dict()
    close_stack = []
    for i in range(n):
      c = db[i]
      if c in self.open_symbols:
        if c not in open_stacks_dict:
          open_stacks_dict[c] = [[]]
        open_stacks = open_stacks_dict[c]
        open_stack = open_stacks[-1]
        if len(open_stack) > 0 and open_stack[-1] + 1 != i:
          open_stacks.append([])
          open_stack = open_stacks[-1]
        open_stack.append(i)
      elif c in self.close_symbols:
        open_symbol = self.get_open_symbol(c)
        if open_symbol not in open_stacks_dict:
          raise Exception("Wrong dot-bracket RNA 2D structure!")

        if len(close_stack) > 0:
          previous_close_index = close_stack[-1]
          previous_close_symbol = db[previous_close_index]
          if previous_close_symbol != c or previous_close_index + 1 != i:
            open_symbol1 = self.get_open_symbol(previous_close_symbol)
            open_stacks1 = open_stacks_dict[open_symbol1]
            open_stack1 = open_stacks1[-1]
            npairs = len(close_stack)
            helices.append(open_stack1[-npairs:] + close_stack)
            if len(open_stack1) == 0:
              if len(open_stacks1) > 1:
                open_stacks1.pop()
              else:
                del open_stacks_dict[open_symbol1]
#                del open_stacks1 # Wrong!!!
            else:
              del open_stack1[-npairs:]
            close_stack = []

        close_stack.append(i)
        open_stacks = open_stacks_dict[open_symbol]
        open_stack = open_stacks[-1]
        if len(open_stack) == len(close_stack):
          helices.append(open_stack + close_stack)
          if len(open_stacks) > 1:
            open_stacks.pop()
          else:
#            del open_stacks  # Wrong!!!
            del open_stacks_dict[open_symbol]
          close_stack = []
    return helices

  def helices_to_graph(self, helices, l):
    G = nx.MultiGraph()
    nodes = set()
    edges = set() # 0: pairing, 1: stacking/helix region, 2: loop
    for h in helices:
      h.sort()
      n = len(h)
      if n == 2:
        nodes.add(h[0])
        nodes.add(h[1])
        edges.add((h[0], h[1], 0))
      else:
        n2 = int(n/2)
        v = [h[0], h[n2-1], h[n2], h[-1]]
        nodes.update(v)
        edges.add((v[0], v[3], 0))
        edges.add((v[0], v[1], 1))
        edges.add((v[1], v[2], 0))
        edges.add((v[2], v[3], 1))
    nodes.add(0)
    nodes.add(l-1)

    nodes = list(nodes)
    nodes.sort()
    n = len(nodes)
    for i in range(n-1):
      a, b = nodes[i:i+2]
      if (a,b,1) in edges:
        continue
      if (a,b,0) in edges and a + 1 == b:
        continue
      edges.add((a,b,2))

#    print('nodes:', nodes)
#    print('edges:', edges)

    for node in nodes:
      G.add_node(node)

    for edge in edges:
      G.add_edge(edge[0], edge[1], type=edge[2])

    return G

  def db_to_graph(self, db):
    helices = self.find_helices(db)
    n = len(db)
    return self.helices_to_graph(helices, n)

  def draw_graph(self, g, outfile):
    pos = nx.spring_layout(g, seed=0)
    nx.draw_networkx_nodes(g, pos)
    nx.draw_networkx_labels(g, pos, labels={n: n for n in g})
    edge_labels = nx.get_edge_attributes(g, "type")
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    for e in g.edges(data='type'):
      if e[2] == 0:
        style = '--'
      else:
        style = '-'
#  nx.draw_networkx_edges(g, pos, edgelist=[edge], style=style)

      ax.annotate("", xy=pos[e[0]], xycoords='data', xytext=pos[e[1]], textcoords='data', arrowprops=dict(arrowstyle="-", ls=style, color="0.5", shrinkA=5, shrinkB=5, patchA=None, patchB=None, connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])),),)

    plt.savefig(outfile)



