.. _bs-form-label:

###############
Bioscript Forms
###############

This section will describe some tricks about **forms** in Bioscript.


*********
FileField
*********

When your plugin needs file(s) as input, you will use a FileField.
Pay attention that **you must** use the class :class:`~tw2.bs.BsFileField`.
It provides for each input file the ability to upload it from the computer (as a *traditional* file input) or to give an URL to the file - which is really useful when files are too big.
If we get back to this example :download:`Second <../../examples/Second.py>` and replace the *TextField* by a *FileField*, we get::

    from bsPlugins import twb
    class MySimpleForm(BaseForm):
        text_file = twb.BsFileField(label="A file : ")
        submit = twf.SubmitButton(id="submit", value="Write it")

***************
Multiple fields
***************

Sometimes you need a parameter that accepts an arbitrary number of values. For instance, a list of files. In order to do that, you must set **multiple** to **the multiple class name** in the parameter description::

     parameters = {'in': [{'id': 'input', 'type': 'file', 'multiple': 'my_inputs'}, ],
              ... }

And then you inherit your *my_inputs* class from :class:`~twb.BsMultiple`.
It will look like this::

    class MySimpleForm(BaseForm):

        param_one = twf.TextField(label="An unique parameter : ")

        class my_inputs(twb.BsMultiple):
            my_files = twb.BsFileField(label='my files', validator=twb.BsFileFieldValidator(required=True))
            my_parameters = twf.TextField(label='my param', validator=twc.Validator(required=True))

Here you can enter one time the parameter *param_one* and multiple times the parameters *my_files* and *my_parameters*. Then in the **__call__** function of the plugin, you can retrieve each parameter like that::

    def __call__(*args, **kw):
        p1 = kw.get('param_one', None)
        my_list_of_files = kw['my_inputs']['my_files']
        my_list_of_parameters = kw['my_inputs']['my_parameters']

        # Here I am sure that 'my_list_of_files' and 'my_list_of_parameters' will contain
        # the same number of parameters and files as I have specified a validator for each of them.



**********
Validation
**********

=================
Simple validation
=================

You can refer to the `toscawidget documentation on validation <http://tw2core.readthedocs.org/en/latest/design/#validation>`_.

======================
BsFileField Validation
======================

As :class:`~twb.BsFielField` is a bit special, is has a special validator, the class :class:`~twb.BsFielFieldValidator`.
You must use it like that::

    from bsPlugins import twb
    class MySimpleForm(BaseForm):
        text_file = twb.BsFileField(label="A file : ", validator=twb.BsFileFieldValidator(required=True))
        submit = twf.SubmitButton(id="submit", value="Write it")

