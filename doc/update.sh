echo "##################
Plugins collection
##################
" > source/content/bsPlugins.rst

grep BasePlugin ../build/lib/bsPlugins/*.py | perl -ne '/\/(\w+)\.py:class (.*)\(Base/; print ".. autoclass:: bsPlugins.$1.$2\n\n"' >> source/content/bsPlugins.rst

make clean html
ssh -t ginger.epfl.ch "sudo /usr/bin/scp -r '"$USER"@"$HOST":"$PWD"/build/html/*' /srv/bsplugins"


