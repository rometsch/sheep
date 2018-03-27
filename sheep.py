
import os
import tempfile
import shutil
import uuid
import tarfile
import xml.etree.ElementTree as ET
import paramset
import parser

def expand_path(path):
    path = os.path.expanduser(path)
    path = os.path.expandvars(path)
    return path

def abs_expand_path(path, base=os.getcwd()):
    path = expand_path(path)
    if path != os.path.abspath(path):
        path = os.path.join(base, path)
    return path

def parse_file_to_lines( filepath ):
    """ Parse a textfile into a list with lines.
    Remove all text following the cahracter #.
    Discard empty lines."""
    lines = []

    with open(filepath, 'r') as f:
        for line in f:
            # Split the line and remove part after the first #
            line = line.strip().split('#')[0].strip()
            if line == "":
                continue
            lines.append(line)
    return lines

def copy(src, dst):
    """ Copy both files and directories from src to dst.
    For files, copy2() and for directories copytree is used. """
    try:
        shutil.copy2(src, dst)
    except IsADirectoryError:
        try:
            shutil.copytree(src, dst)
        except FileExistsError:
            shutil.copytree(src, os.path.join(dst, os.path.basename(src)))

class Sheep:

    def __init__(self, setup_dir, config_file='sheep.xml'):
        self.uuid = str(uuid.uuid4())
        self.temp_dir = None
        self.temp_dir_obj = None
        self.src_list = []
        self.scripts = {}
        self.param_set = None
        self.parameters = {}
        self.parameters_in_config = {}
        self.tar_file = None
        self.setup_dir = os.path.abspath(setup_dir)
        # If the config dir is not an absolute path assume its meant to be inside the setup dir.
        if config_file != os.path.abspath(config_file):
            config_file = os.path.join( setup_dir, config_file )
        self.cfg = ET.parse(config_file);
        self.parse_src_list()
        #self.parse_parameter_config()
        self.provide_temp_dir()
        self.copy_src()
        self.copy_scripts()
        self.construct_param_set()


    def provide_temp_dir(self):
        """ Create a save temporary directory, unique for every instance of sheep. """
        if self.temp_dir is None:
            self.temp_dir_obj = tempfile.TemporaryDirectory(prefix = self.uuid)
            self.temp_dir = self.temp_dir_obj.name

    def parse_src_list(self):
        """ Provide a list of paths to be copied to the temp directory
        directory. User variables and the ~ shorthand are expanded. """
        self.src_list = [c.text for c in self.cfg.find('./source')]

    def copy_scripts(self):
        """ Provide a list of paths of scripts. They need to be copied aswell.
        User variables and the ~ shorthand are expanded. """
        scripts = self.cfg.find('./scripts')
        for script in scripts:
            path = script.text
            self.scripts[script.tag] = path
            copy( abs_expand_path(path, base = self.setup_dir),
            os.path.join(self.temp_dir, "{}_sheep".format(script.tag) ) )

    def parse_parameter_config(self, parameters_file):
        """ Load the names of the parameters, how to translate them from generic names
        and the desitination file from the parameters.txt file. """
        lines = parse_file_to_lines( parameters_file )
        for line in lines:
            parts = line.split()
            if parts[1] == '-':
                parts[1] = parts[0]
            self.parameters[parts[0]] = parts[1:]

    def get_temp_path(self, filename ):
        """ Return the absolute path of the file *filename* inside the temp dir. """
        return os.path.join( self.temp_dir, filename)

    def list_temp_dir(self):
        """ Return a list with absolute paths of every file or directory in
        the temp folder. """
        return [os.path.join(self.temp_dir, f) for f in os.listdir(self.temp_dir)]

    def copy_src(self):
        if len(self.src_list) == 0:
            print("Warning: Nothing is copyied b.c. src_list is empty.")
        for path in self.src_list:
            try:
                copy(abs_expand_path(path, base = self.setup_dir), self.temp_dir)
            except TypeError:
                print("Error while trying to copy {}".format(path));


    def make_tar(self):
        """ Make a gzipped tar file containing all files and directories
        inside the tmp folder. """
        self.save_changes()
        self.tar_file = self.get_temp_path( 'content.tgz' )
        files = self.list_temp_dir()
        with tarfile.open( self.tar_file , 'x:gz') as tf:
            for path in files:
                tf.add(path, arcname = path.replace(self.temp_dir,'').lstrip('/'))

    def save_tar(self, dst):
        """ Construct and move the tar file to the given location. """
        if self.tar_file is None:
            self.make_tar()
        shutil.move(self.tar_file, dst)

    def enforce_param_known(self, name):
        if not (name in self.parameters or name in self.parameters_in_config):
            raise ValueError("Setting the parameter {} is not allowed. Make sure it is spelled correct and that its either definde in the config or already present inside a config file.".format(name))

    def construct_param_set(self):
        """ Make a param set object using the parser type specified in the config """
        parserType = self.cfg.find('./paramset/type').text
        paramFile = self.cfg.find('./paramset/path').text
        paramFile = self.get_temp_path(paramFile)
        try:
            Parser = parser.avail[parserType]
            par = Parser(paramFile)
            self.param_set = paramset.ParamSet(par)
        except KeyError:
            print("Could not find parser type {}".format(parserType))

    def set_param(self, param, value ):
        try:
            self.param_set.set_param(param, value)
        except KeyError:
            print("Could not find parameter {}".format(param))
            raise

    def save_changes(self):
        """ Write saves to disk """
        if not self.param_set is None:
            self.param_set.parser.save()

    def translate_param_name(self, to_convention, name ):
        pass
