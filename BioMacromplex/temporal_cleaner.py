import os
import shutil
if not os.path.exists('tmp'):
    os.mkdir('tmp')


def clean_temporal():
    shutil.rmtree('tmp/')

