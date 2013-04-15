.. _bs-form-label:

###############
Bioscript Forms
###############

This section will describe some tricks about **forms** in Bioscript.


*********
FileField
*********

When your plugin needs file(s) as input, you will use FileField.
Pay attention that **you must** use the class :class:`~tw2.bs.BsFileField`.
They provide for each users the ability to, when a file is asked as input either provide a file form
their own computer or provide an URL which link to the file. Which is really useful when files are too big.
If we get back on our example :download:`file <../../examples/Second.py>` and replace the *TextField* by a *FileField*::


    from bsPlugins import twb
    class MySimpleForm(BaseForm):
        text_file = twb.BsFileField(label="A file : ")
        submit = twf.SubmitButton(id="submit", value="Write it")

***************
Multiple fields
***************

Sometimes you need to input a parameter multiple times. For instance a list of files and parameters. In order to do that, you must
set **multiple true** in the parameter description::

     parameters = {'in': [{'id': 'input', 'type': 'file', 'multiple': True}, ],
              ... }

And then you must use the class :class:`~twb.BsMultiple`.
It will look like this::

    class MySimpleForm(BaseForm):

        param_one = twf.TextField(label="An unique parameter : ")

        class my_inputs(twb.BsMultiple):
            my_files = twb.BsFileField(label='my files', validator=twb.BsFileFieldValidator(required=True))
            my_parameters = twf.TextField(label='my param', validator=twc.Validator(required=True))

Here you can enter one time the parameter *param_one* and multiple times the parameters *my_files* and *my_parameters*.
The in the **__call__** function of the plugin, you can retrieve each parameters like that::

    def __call__(*args, **kw):
        p1 = kw.get('param_one', None)
        my_list_of_files = kw['my_inputs']['my_files']
        my_list_of_parameters = kw['my_inputs']['my_parameters']

        # Here I'm sure that 'my_list_of_files' & 'my_list_of_parameters' will contains
        # the same amount of parameters and files as I have specified a validator
        # on each of them.



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

As :class:`~twb.BsFielField` are a bit specials, they have a special validator, the class :class:`~twb.BsFielFieldValidator`.
You must use like that::

    from bsPlugins import twb
    class MySimpleForm(BaseForm):
        text_file = twb.BsFileField(label="A file : ", validator=twb.BsFileFieldValidator(required=True))
        submit = twf.SubmitButton(id="submit", value="Write it")



