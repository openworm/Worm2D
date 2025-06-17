/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
/*
 * owINeuronSimulator.h
 * Interface for simulator of neuronal activity
 *  Created on: Jun 27, 2016
 *      Author: Sergey Khayrulin
 */

#include <vector>
#include <Python.h>

#ifndef INC_OWINEURONSIMULATOR_H_
#define INC_OWINEURONSIMULATOR_H_



class owINeuronSimulator {
protected:
	std::vector<float> unpackPythonList(PyObject* pValue, size_t musclesNum=96){
		Py_ssize_t size = PyList_Size(pValue);
		//std::vector<float> test(musclesNum); //needs to change! 96 is hardcoded
		std::vector<float> test((size_t) size);
		//printf("====\n");
		for (Py_ssize_t i = 0; i < size; i++) {
			float value;
			value = (float)PyFloat_AsDouble(PyList_GetItem(pValue, i));
			test[i]= value;
		}
		Py_DECREF(pValue);
		return test;
	}
	PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *pClass, *pInstance, * nrn_sim;
public:
	//virtual std::vector<float> run() = 0;
	virtual void run() = 0;
	virtual ~owINeuronSimulator()
	{
		Py_DECREF(pName);
		Py_DECREF(pModule);
		Py_DECREF(pDict);
		Py_DECREF(pFunc);
		Py_DECREF(pValue);
		Py_DECREF(pClass);
		Py_DECREF(pInstance);
		Py_DECREF(nrn_sim);
		Py_FinalizeEx();
	}
};

#endif /* INC_OWINEURONSIMULATOR_H_ */
