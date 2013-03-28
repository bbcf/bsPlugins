###########
# IMPORTS #
###########
import os
import random
from bsPlugins import OperationPlugin
from bsPlugins import BaseForm
from bsPlugins import twf

############
# METADATA #
############

meta = {'version': "1.0.1",
        'author': "Yohan Jarosz",
        'contact': "webmaster-bbcf@epfl.ch"}

########
# FORM #
########


class MySimpleForm(BaseForm):
    text = twf.TextField(label="Input something in a file : ")
    submit = twf.SubmitButton(id="submit", value="Write it")

######################
# PLUGIN INFORMATION #
######################
parameters = {'in': [{'id': 'input', 'type': 'text', 'required': True}, ],
              'out': [{'id': 'output', 'type': 'file'}, ]}

plugin_information = {
    'title': 'WriteFile',

    'description': """As an exemple, this plugin writes the input you give in a text file.
    You can also put <a href='https://github.com/bbcf/bs.operations/blob/master/Examples.py'>links</a> to some documentation or
    other <br/> HTML tags.""",

    'path': ['Examples', 'Automatics forms', 'Write output to a file'],

    'meta': meta,

    'in': parameters['in'],
    'out': parameters['out'],

    'output': MySimpleForm
}



##########
# PLUGIN #
##########


class Simple(OperationPlugin):

    info = plugin_information

    def __call__(self, *args, **kw):
        text = kw.get('input', '')                    # get the parameter back

        path = self.temporary_path()                  # get a temporary path
        with open(path, 'w') as f:                    # open a file & write the input
            f.write(text)
        self.new_file(path, 'output')                 # add a file to the result

        return self.display_time()                    # display the time elapsed











from bsPlugins import OperationPlugin, BaseForm, DynForm
from bsPlugins import twf, twc, twd
# common information for the two plugins
meta = {'version': "1.0.0",
        'author': "Yohan Jarosz",
        'contact': "webmaster-bbcf@epfl.ch"}

###################################################### PLUGIN ONE
# import toscawidget2 modules in order to build forms


class OutputForm(BaseForm):
    # the parameter 'input'
    input = twf.FileField(label_text="My input", validator=twf.FileValidator(required=True))

    # the submit button
    submit = twf.SubmitButton(id="submit", value="Submit My job")


class WithForm(OperationPlugin):

    info = {
        'title': 'Simple Form customization',
        'description': """See <a href="http://tw2core.readthedocs.org/en/latest/index.html">
        toscawidget documentation</a> for more information.""",
        'path': ['Examples', 'User defined forms', 'Static'],
        'in': [{'id': 'input', 'type': 'file', 'required': True}],
        'out': [],
        'meta': meta,
        'output': OutputForm                         # Define the form you want to use
    }

    def __call__(self, *args, **kw):     # proceed as usual
        file_input = kw.get('input', '')
        return file_input

###################################################### PLUGIN TWO


class OutputForm(BaseForm):
    # the parameter 'input'
    input = twf.FileField(label_text="My input not required")

    # the submit button
    submit = twf.SubmitButton(id="submit", value="Submit My job")


class WithAnotherForm(OperationPlugin):

    info = {
        'title': 'Another form exemple',
        'description': """See <a href="http://tw2core.readthedocs.org/en/latest/index.html">
        toscawidget documentation</a> for more information.""",
        'path': ['Examples', 'User defined forms', 'Static 2'],
        'in': [{'id': 'input', 'type': 'file', 'required': True}],
        'out': [],
        'meta': meta,
        'output': OutputForm                         # Define the form you want to use
    }

    def __call__(self, *args, **kw):     # proceed as usual
        file_input = kw.get('input', '')
        return file_input

###################################################### PLUGIN THREE
# import dynamic modules


# class DynamicOutputForm(DynForm):
#     # wrap dynamic content with a HidindTableLayout for exemple
#     class method(twd.HidingTableLayout):
#         method = twd.HidingSingleSelectField(label='Select method', options=('This is', 'a-demo'),
#             mapping={
#                 'This is': ['one'],
#                 'a-demo': ['two', 'three'],
#                     })
#         one = twf.TextField(label='One')
#         two = twf.TextField(label='Two')
#         three = twf.TextField(label='Three')

#     class parameters(twd.HidingTableLayout):
#         params = twd.HidingSingleSelectField(label='Select parameters', options=('opt1', 'opt2', 'opt3'),
#             mapping={
#                 'opt1': ['a'],
#                 'opt2': ['a', 'b'],
#                 'opt3': ['a', 'b', 'c'],
#                     })
#         a = twf.TextField(label='1')
#         b = twf.TextField(label='2')
#         c = twf.TextField(label='3')

#     #the submit button
#     submit = twf.SubmitButton(id="submit", value="Submit My job")


# class WithDynamicForm(OperationPlugin):

#     info = {
#         'title': 'Dynamic form',
#         'description': """See <a href="http://tw2core.readthedocs.org/en/latest/index.html">
#         toscawidget documentation</a> for more information.""",
#         'path': ['Examples', 'User defined forms', 'Dynamic'],
#         'in': [{'id': 'input', 'type': 'text', 'required': True}],
#         'out': [],
#         'meta': meta,
#         'output': DynamicOutputForm                         # Define the form you want to use
#     }

#     def __call__(self, *args, **kw):     # proceed as usual
#         print 'got args %s & kw %s' % (args, kw)
#         # get method back :
#         method_selected = kw.get('method:method')
#         method_one = kw.get('method:one')
#         method_two = kw.get('method:one')
#         method_three = kw.get('method:one')
#         # get parameters back
#         params_selected = kw.get('parameters:param')
#         a = kw.get('parameters:a')
#         b = kw.get('parameters:b')
#         c = kw.get('parameters:c')

#         return self.display_time()
