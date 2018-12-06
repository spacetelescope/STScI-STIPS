"""
General code for handling the segmented data tables that STIPS uses in order to avoid issues very
long tables.

:Author: Brian York

:Organization: Space Telescope Science Institute

"""
from __future__ import absolute_import, division, print_function

# External modules
import os, sys
from astropy.io import ascii, fits
from astropy.table import Table

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from cStringIO import StringIO

class StipsDataTable(object):
    def __init__(self, **kwargs):
        """
        This function creates a data table object. It generally expects to have a file name, and 
        will take formatting arguments (and the like), and do both input and output.
        """
        self._file = kwargs.get('file_name', None)
        self._in_dir = kwargs.get('in_dir', os.getcwd())
        self._format = kwargs.get('format', 'ascii.ipac')
        self._chunk_size = kwargs.get('chunk_size', 100000)
        self.init_table()
    
    @classmethod
    def initFromFile(cls, table_file, **kwargs):
        """
        Creates a data table given a file
        """
        (file_path, file_name) = os.path.split(table_file)
        (table_base, ext) = os.path.splitext(file_name)
        if ext == 'fits':
            format = 'fits'
        else:
            format = 'ascii.ipac'
        table = cls(file_name=table_file, format=format, in_dir=file_path, **kwargs)
        return table
    
    @property
    def file(self):
        return self._file
    
    @property
    def in_dir(self):
        return self._in_dir
    
    @property
    def format(self):
        return self.format
    
    @property
    def chunk_size(self):
        return self._chunk_size
    
    def init_table(self):
        """
        Nothing to do in base class
        """
        self._current_chunk = 0
    
    def read_chunk(self, chunk=None):
        raise NotImplementedError("Read Chunk not implemented in base class")
    
class StipsFitsTable(StipsDataTable):
    def init_table(self):
        super(StipsFitsTable, self).init_table()
        if os.path.isfile(self.file):
            t = self.read_chunk()
            self.meta = t.meta
    
    def read_chunk(self, chunk=None):
        if chunk is None:
            chunk = self._current_chunk
        return Table.read(self.file, hdu=chunk)
    
    def write_chunk(self, chunk):
        if os.path.isfile(self.file):
            hdulist = fits.HDUList.fromfile(self.file)
        else:
            hdulist = fits.HDUList()
        new_hdu = chunk.table_to_hdu()
        hdulist.insert(self._current_chunk, new_hdu)
        self._current_chunk += 1
        hdulist.writeto(self.file, overwrite=True)
        
class StipsIpacTable(StipsDataTable):
    def init_table(self):
        super(StipsIpacTable, self).init_table()
        self.header_length = 4
        if os.path.isfile(self.file):
            self.names, self.meta = self.init_names()
    
    def read_chunk(self, chunk=None)
        if chunk is None:
            chunk = self._current_chunk
        lines = []
        with open(self.file) as inf:
            for i in range(self.chunk_size*chunk):
                inf.readline()
            for i in range(self.chunk_size):
                lines.append(inf.readline().rstrip(os.linesep))
        return Table.read(lines, format='no_header', names=self.names, guess=False)
    
    def write_chunk(self, chunk):
        if not os.path.isfile(self.file):
            with open(self.file, 'w') as outf:
                outf.write("\n")
            data = StringIO()
            chunk.write(data, format='ascii.ipac')
            data.seek(0)
            if self._current_chunk != 0:
                for i in range(self.header_length):
                    data.readline()
            self._current_chunk += 1
            with open(self.file, 'a') as outf:
                outf.write(data.read())
        
    def init_names(self):
        """
        Initialize names for the ascii.ipac table
        """
        names = None
        lines = []
        with open(self.file) as inf:
            for i, line in enumerate(inf):
                lines.append(line.rstrip(os.linesep))
                if i == self.chunk_size:
                    break
        t = Table.read(lines, format=self.format, guess=False)
        return (chunk.colnames, chunk.meta)

   
