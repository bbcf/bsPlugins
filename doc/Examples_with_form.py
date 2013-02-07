from bs.operations import base  # import the base class to build your plugin

# common information for the two plugins
meta = {'version': "1.0.0",
        'author': "Yohan Jarosz",
        'contact': "webmaster-bbcf@epfl.ch"}

###################################################### PLUGIN ONE
# import toscawidget2 modules in order to build forms
import tw2.forms as twf


class OutputForm(base.BaseForm):
    # the parameter 'input'
    input = twf.FileField(label_text="My input", validator=twf.FileValidator(required=True))

    # the submit button
    submit = twf.SubmitButton(id="submit", value="Submit My job")


class WithForm(base.OperationPlugin):

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


class OutputForm(base.BaseForm):
    # the parameter 'input'
    input = twf.FileField(label_text="My input not required")

    # the submit button
    submit = twf.SubmitButton(id="submit", value="Submit My job")


class WithAnotherForm(base.OperationPlugin):

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
import tw2.dynforms as twd  # import dynamic modules


class DynamicOutputForm(base.DynForm):
    # wrap dynamic content with a HidindTableLayout for exemple
    class method(twd.HidingTableLayout):
        method = twd.HidingSingleSelectField(label='Select method', options=('This is', 'a-demo'),
            mapping={
                'This is': ['one'],
                'a-demo': ['two', 'three'],
                    })
        one = twf.TextField(label='One')
        two = twf.TextField(label='Two')
        three = twf.TextField(label='Three')

    class parameters(twd.HidingTableLayout):
        params = twd.HidingSingleSelectField(label='Select parameters', options=('opt1', 'opt2', 'opt3'),
            mapping={
                'opt1': ['a'],
                'opt2': ['a', 'b'],
                'opt3': ['a', 'b', 'c'],
                    })
        a = twf.TextField(label='1')
        b = twf.TextField(label='2')
        c = twf.TextField(label='3')

    #the submit button
    submit = twf.SubmitButton(id="submit", value="Submit My job")


class WithDynamicForm(base.OperationPlugin):

    info = {
        'title': 'Dynamic form',
        'description': """See <a href="http://tw2core.readthedocs.org/en/latest/index.html">
        toscawidget documentation</a> for more information.""",
        'path': ['Examples', 'User defined forms', 'Dynamic'],
        'in': [{'id': 'input', 'type': 'text', 'required': True}],
        'out': [],
        'meta': meta,
        'output': DynamicOutputForm                         # Define the form you want to use
    }

    def __call__(self, *args, **kw):     # proceed as usual
        print 'got args %s & kw %s' % (args, kw)
        # get method back :
        method_selected = kw.get('method:method')
        method_one = kw.get('method:one')
        method_two = kw.get('method:one')
        method_three = kw.get('method:one')
        # get parameters back
        params_selected = kw.get('parameters:param')
        a = kw.get('parameters:a')
        b = kw.get('parameters:b')
        c = kw.get('parameters:c')

        return self.display_time()
