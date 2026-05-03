import os
import ttbarEFT

pjoin = os.path.join

# This function takes as input any path (inside of ttbarEFT/ttbarEFT), and returns the absolute path
def ttbarEFT_path(path_in_repo):
    return pjoin(ttbarEFT.__path__[0], path_in_repo)
