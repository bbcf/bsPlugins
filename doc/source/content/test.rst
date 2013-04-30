#############
Test a plugin
#############

 Yes, please test them... And make the test quickly reproducible (search for "unit tests").

*************
Test the form
*************

Navigate to `bsPLugins/tests <https://github.com/bbcf/bsPlugins/tree/master/tests>`_ and in the terminal type (replace *<PluginName>* by the name of your plugin)::

 python test_form.py

It should output::

 Test plugin <PluginName>
 serving on http://127.0.0.1:8080

Then open your browser at the latter address. It should display the form as it will be for Bioscript users. After you enter your data, click on the submit button: it will print a JSON (a string representing a dictionary) which is the data as it will be passed to the plugin part (remember the __call__() method ?).

If there are any errors, take time to debug them with the precious help of the marvellous ToscaWidgets documentation.

******************
Test the operation
******************

=======
By Hand
=======

Although it is a class, it acts as a function since it has a **__call__** method. So you can call it as MypluginPlugin()(<arguments>) - note the pair of parenthesis before the arguments, to create an instance of the class. For instance::

 SmoothingPlugin()(track=..., assembly='mm9', ...)

I advise to instanciate separately though, in order to reuse the plugin instance later::

 plugin = SmoothingPlugin()               # instanciate
 plugin(track=..., assembly='mm9', ...)   # start the operation

However, the best way to test it is to use a dictionary (such as the JSON the form outputs), this way::

 kw = {'track':..., 'assembly':'mm9', ...}
 plugin(**kw)

This is what you will see in unit tests that already exist (and that you are encouraged to imitate) in `bsPLugins/tests <https://github.com/bbcf/bsPlugins/tree/master/tests>`_ , e.g. test_Smoothing.py.

Ouput files are retrived as follows::

 plugin.output_files  # a list of lists: [[filename,type], [filename,type], ...]

So to get the first (and often unique) output file path just do::

 plugin.output_files[0][0]

===========
PRINT DEBUG
===========

When a plugin is running on the web interface, if an operation fails, you cannot really see the print statements you've set, but you have a method to do that::
    
 self.debug(<some debug statement>)

You will see the output in the 'full traceback' section. This is available on the <dev> server.


..note::
 
 self.debug() will also print the <debug statement> to the standard out if you set self.is_debug = True.

============
Unit testing
============

You can run already existing unit tests by going to this folder and typing in the console::

 nosetests test_Smoothing.py

If it does not know **nosetests** ou **unittest2**, install it first::

 pip install nose unittest2

To disable ToscaWidgets debug warnings, you may prefer to call::

 nosetests --logging-filter=-tw2 test_Smoothing.py

And make your own !!

