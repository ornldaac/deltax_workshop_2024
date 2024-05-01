#! /usr/bin/env python
from __future__ import print_function
import sys
import shutil
import os

# Install attached files in appropriate locations
# Find anuga:
try:
    import anuga
except:
    print('\n')
    print('#' * 60)
    print('CANNOT INSTALL THE FILES')
    print('ANUGA must be installed before running this script.')
    print('Find it at https://github.com/GeoscienceAustralia/anuga_core')
    print('#' * 60)
    print('\n')
    sys.exit()
    
path_to_anuga = anuga.__path__[0]

# Files to copy
file1 = 'baptist_operator.py'
file2 = 'config.py'
file3 = 'shallow_water_domain.py'

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)

# Copy files
src1 = os.path.join(path,file1)
dst1 = os.path.join(path_to_anuga,'operators',file1)
shutil.copy2(src1, dst1)

src2 = os.path.join(path,file2)
dst2 = os.path.join(path_to_anuga,file2)
#shutil.copy2(src2, dst2)

src3 = os.path.join(path,file3)
dst3 = os.path.join(path_to_anuga,'shallow_water',file1)
shutil.copy2(src3, dst3)

# Success message
print('Success!')
print('Installed extra codes in appropriate folders')
