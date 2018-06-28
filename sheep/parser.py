import os
from collections import OrderedDict as ODict

def parse_to_dict(lines):
    # Initialize a first group to store key value pairs listed before the first group.
    groups = ODict()
    stats = ODict()
    name = '_root'
    groups[name] = ODict()
    stats[name] = [0,0]
    num_comment = 1
    for line in lines:
        line = line.strip()
        # Skip empty lines
        if line == "":
            continue
        # Identify a header line
        if line[0] == '[':
            if line[-1] != ']':
                raise ValueError('Found a line starting wih "[" but not ending in "]"!')
            name = line.lstrip('[').rstrip(']')
            groups[name] = ODict()
            stats[name] = [0, 0]
        # Append normal lines to the last group
        else:
            if line[0] in ["#",";"]:
                groups[name]["__comment_" + "{}".format(num_comment)] = line
                num_comment += 1
            else:
                parts = line.split()
                groups[name][parts[0]] = parts[1:]
                stats[name] += [0]*(len(parts) - len(stats[name]))
                for n,p in enumerate(parts):
                    stats[name][n] = max(stats[name][n], len(p))


    return (groups, stats)

def parse_dict_to_str(dct, stats=None):
    lines = []
    for group in dct:
        # blank line before the group header
        if len(lines) > 0:
            lines.append('')
        # Write the header if group is not the root.
        if group != '_root':
            lines.append('[{}]'.format(group))
            lines.append('')
        # Write the key value pairs
        for key in dct[group]:
            if key[:9] == "__comment":
                lines.append("")
                lines.append(dct[group][key])
                lines.append("")
            elif stats is None:
                lines.append('{}   {}'.format( key, '   '.join(dct[group][key])  ))
            else:
                st = stats[group]
                line = format(key, '{}s'.format(st[0]))
                for n, val in enumerate(dct[group][key]):
                    line += '   ' + format(val, '{}s'.format(st[n+1]))
                lines.append(line)
    return lines

class IniParser:

    def __init__(self, inifile):
        self.inifile = os.path.abspath(inifile)
        self.load()

    def load(self):
        """ Load the ini file into memory. """
        with open(self.inifile, 'r') as rf:
            self.cfg_dct, self.stats = parse_to_dict(rf)

    def parse_dict_to_str(self):
        """ Write the dict into lines. """
        return parse_dict_to_str(self.cfg_dct, self.stats)

    def save(self, path=None):
        """ Write the config back to the file by default or to any file specified by *path*. """
        if path is None:
            path = self.inifile
        with open(path, 'w') as of:
            of.write('\n'.join(self.parse_dict_to_str()))

    def get_param_dict(self):
        """ Return the parameter dict. """
        return self.cfg_dct

# List parser classes available in this implementation.
# This is handy for specifying a parser class in a config file.
avail = {
    'plutoIni' : IniParser,
    'fargo3dIni' : IniParser,
    'fargoIni' : IniParser
}
