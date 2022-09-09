"""
General code for handling the segmented data tables that STIPS uses in order to avoid issues very
long tables.

:Author: Brian York

:Organization: Space Telescope Science Institute

"""
from __future__ import absolute_import, division, print_function

# External modules
import os
from astropy.io import ascii, fits
from astropy.table import Table
from collections import OrderedDict
from io import StringIO


class StipsDataTable(object):
    def __init__(self, **kwargs):
        """
        This function creates a data table object. It generally expects to have a file name, and
        will take formatting arguments (and the like), and do both input and output.
        """
        self._file = kwargs.get('file_name', None)
        self._fname = os.path.split(self._file)[1]
        self._in_dir = kwargs.get('in_dir', os.getcwd())
        self._format = kwargs.get('format', 'ascii.ipac')
        self._chunk_size = kwargs.get('chunk_size', 100000)
        starting_chunk = kwargs.get('starting_chunk', 0)
        self.init_table(starting_chunk)

    @staticmethod
    def dataTableFromFile(table_file, **kwargs):
        """
        Creates a data table given a file
        """
        table_classes = {'fits': StipsFitsTable, 'ipac': StipsIpacTable, 'default': StipsIpacTable}
        (file_path, file_name) = os.path.split(table_file)
        (table_base, ext) = os.path.splitext(file_name)
        if 'fits' in ext:
            table = table_classes['fits'](file_name=table_file, format='fits', in_dir=file_path, **kwargs)
        else:
            table = table_classes['default'](file_name=table_file, format='ascii.ipac', in_dir=file_path, **kwargs)
        return table

    @property
    def file(self):
        return self._file

    @property
    def in_dir(self):
        return self._in_dir

    @property
    def format(self):
        return self._format

    @property
    def chunk_size(self):
        return self._chunk_size

    @property
    def chunk(self):
        return self._current_chunk

    def init_table(self, starting_chunk):
        """
        Nothing to do in base class
        """
        self._current_chunk = starting_chunk

    def read_chunk(self, chunk=None, advance=True, set=True):
        raise NotImplementedError("Read Chunk not implemented in base class")

    def write_chunk(self, chunk=None, advance=True):
        raise NotImplementedError("Write Chunk not implemented in base class")


class StipsFitsTable(StipsDataTable):
    def init_table(self, starting_chunk):
        super(StipsFitsTable, self).init_table(starting_chunk)
        if os.path.isfile(self.file):
            t = self.read_chunk(chunk=0, set=False, advance=False)
            self.meta = t.meta
            self.columns = t.columns

    def read_chunk(self, chunk=None, advance=True, set=True):
        if chunk is None:
            chunk = self._current_chunk
        with fits.open(self.file, memmap=True) as ff:
            if chunk >= len(ff):
                return None
            while not isinstance(ff[chunk], fits.BinTableHDU):
                chunk += 1
        if set:
            self._current_chunk = chunk
        table_data = Table.read(self.file, hdu=chunk)
        if advance:
            self._current_chunk += 1
        return table_data

    def write_chunk(self, chunk, advance=True):
        if chunk is not None:
            if os.path.isfile(self.file):
                hdulist = fits.HDUList.fromfile(self.file)
            else:
                hdulist = fits.HDUList()
            if chunk.meta is None or len(chunk.meta) == 0:
                chunk.meta = self.meta
            try:
                new_hdu = fits.table_to_hdu(chunk)
            except Exception as e:
                raise e
            hdulist.insert(self._current_chunk, new_hdu)
            if advance:
                self._current_chunk += 1
            hdulist.writeto(self.file, overwrite=True)


class StipsIpacTable(StipsDataTable):
    def init_table(self, starting_chunk):
        super(StipsIpacTable, self).init_table(starting_chunk)
        self.header_length = 4
        if os.path.isfile(self.file):
            self.columns, self.meta = self.read_metadata(self.file)

    @staticmethod
    def read_metadata(filename, n_lines=100, format='ascii.ipac'):
        lines = []
        for i, line in enumerate(open(filename, 'r')):
            lines.append(line.rstrip(os.linesep))
            if i == n_lines:
                break
        t = Table.read(lines, format=format)
        if "keywords" in t.meta:
            meta = OrderedDict()
            for item in t.meta["keywords"]:
                meta[item] = t.meta["keywords"][item]["value"]
            t.meta = meta
        return t.columns, t.meta

    @staticmethod
    def read_table(filename, n_chunk=100000, format="ipac"):
        """
        Chunk reader to (hopefully) not take up so much memory when reading very large tables.
        """
        names = None
        lines = []
        for i, line in enumerate(open(filename, 'r')):
            lines.append(line.rstrip(os.linesep))
            if i % n_chunk == n_chunk - 1:
                if i < n_chunk:  # first
                    chunk = ascii.read(lines, format=format, guess=False)
                    names, meta = StipsIpacTable.read_metadata(filename, format=format)
                    chunk.meta = meta
                    yield chunk
                else:
                    chunk = ascii.read(lines, format='no_header', names=names, guess=False)
                    names, meta = StipsIpacTable.read_metadata(filename, format=format)
                    chunk.meta = meta
                    yield chunk
                lines = []
        if lines:
            if names is not None:
                chunk = ascii.read(lines, format='no_header', names=names, guess=False)
                names, meta = StipsIpacTable.read_metadata(filename, format=format)
                chunk.meta = meta
                yield chunk
            else:
                chunk = ascii.read(lines, format=format, guess=False)
                names, meta = StipsIpacTable.read_metadata(filename, format=format)
                chunk.meta = meta
                yield chunk

    def read_chunk(self, chunk=None, advance=True, set=True):
        if not hasattr(self, "_iterator"):
            self._iterator = self.read_table(self.file, self.chunk_size)
        return next(self._iterator, None)

    def write_chunk(self, chunk, advance=True):
        if chunk is not None:
            if not os.path.isfile(self.file):
                with open(self.file, 'w') as outf:
                    outf.write("\\ {}\n".format(self.meta.get('name', 'Default Table')))
                    outf.write("\\ \n")
                    outf.write("\\ Parameters:\n")
                    outf.write("\\ \n")
                    for key in self.meta:
                        if key != 'name':
                            if isinstance(self.meta[key], str):
                                outf.write('\\{}="{}"\n'.format(key, self.meta[key]))
                            else:
                                outf.write('\\{}={}\n'.format(key, self.meta[key]))
                    outf.write("\\ \n")
            data = StringIO()
            if chunk.meta is None or len(chunk.meta) == 0:
                chunk.meta = self.meta
            chunk.write(data, format='ascii.ipac')
            data.seek(0)
            if self._current_chunk != 0:
                for i in range(self.header_length):
                    data.readline()
            if advance:
                self._current_chunk += 1
            with open(self.file, 'a') as outf:
                outf.write(data.read())

    def init_names(self):
        """
        Initialize names for the ascii.ipac table
        """
        lines = []
        with open(self.file) as inf:
            for i, line in enumerate(inf):
                lines.append(line.rstrip(os.linesep))
                if i == self.chunk_size:
                    break
        return self.read_metadata(self.file, format=self.format)
