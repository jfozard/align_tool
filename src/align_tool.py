


"""

Tool for aligning multiple confocal stacks (tiles, timepoints, staged samples) in 2.5 D

Model

"""

from __future__ import print_function


import numpy as np
import numpy.linalg as la
import sys
import math
import itertools
import heapq
import os
import argparse

import time

import scipy.ndimage as nd

from import_tiff import load_tiff

from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMessageBox
            

def scale(x):
    return np.diag(np.append(x,1))

def to_qimg(m):
    m = np.ascontiguousarray(np.dstack([m,m,m]))
    height, width, channels = m.shape
    bytesPerLine = width*3
    return QtGui.QImage(QtCore.QByteArray(m.tostring()), width, height, bytesPerLine, QtGui.QImage.Format_RGB888)
                                                                                                            


class CoordinateTransform(object):
    name = "CoordinateTransform"

    def get_name(self):
        return cls.name
    


class AffineTransform(CoordinateTransform):
    name = "AffineTransform"

    def __init__(self, matrix):
        self.matrix = matrix

    def get_inverse(self):
        return AffineCoordinateTransform(la.inv(self.matrix))

    def map_points(self, x):
        return np.dot(self.matrix[:3,:3], x) + self.matrix[:,:3][np.newaxis,:]

    


"""
Confocal stack - for the moment always acquired in a regular cuboidal coordinate system
"""
class Stack(object):
    def __init__(self, data, spacing):
        self.data = data
        self.spacing = spacing
        self.shape = data.shape
        self.default_transform = AffineTransform(scale(self.spacing))
        self.user_transform = AffineTransform(np.eye(4))
        self.make_maxproj()
        self.visible = True
        self.transparent = False
        
    def make_maxproj(self):
        self.maxproj = np.max(self.data, axis=0)


    def get_transform2d(self):
        m = self.user_transform.matrix.dot(self.default_transform.matrix)
        return m[1:,1:].T.flatten()
    
class WorldController(object):
    def __init__(self):
        self.stacks = []
        self.sbar = None
        self.rw = None
        
    def load_stack(self, filename):
        data, spacing = load_tiff(filename)
        print('load_stack', filename, data.shape, spacing)
        self.stacks.append(Stack(data, spacing[::-1]))
        if self.sbar:
            self.sbar.addObject(filename, self.stacks[-1])
        if self.rw:
            self.rw.add_stack(self.stacks[-1])
       

class WorldViewController(object):
    def __init__(self):
        pass

class StackSideBar(QtWidgets.QDockWidget):
    def __init__(self):
        QtWidgets.QDockWidget.__init__(self)
        splitter = QtWidgets.QSplitter(self)
        splitter.setOrientation(QtCore.Qt.Vertical)
        self.setWidget(splitter)

        #        self.list_widget = QtGui.QTreeWidget()
        self.list_widget = QtWidgets.QTreeWidget()
        self.list_widget.setColumnCount(1)
        self.list_widget.setHeaderLabels(["Name"])
        
        splitter.addWidget(self.list_widget)
        self.list_widget.itemChanged.connect(self.toggleVis)
        self.list_widget.itemSelectionChanged.connect(self.itemSelectionChanged)
    
        self.label = QtWidgets.QLabel("Object details")
        splitter.addWidget(self.label)


    def addObject(self, name, obj):
#        item = QtGui.QListWidgetItem(filename)
        item = QtWidgets.QTreeWidgetItem(None)
        item.setText(0, name)
        item.obj = obj
        item.obj_name = name
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(0, QtCore.Qt.Checked)
        item.prevstate = item.checkState(0)

        if True: #obj.parent == None:
            self.list_widget.addTopLevelItem(item)
        else:
            l = self.getObjItem(obj.parent)
            if l:
                l.addChild(item)
            else:
                self.list_widget.addTopLevelItem(item)

    def removeObject(self, obj):
        for i in range(self.list_widget.topLevelItemCount()):
            l = self.list_widget.topLevelItem(i)
            if l.obj == obj:
                self.list_widget.takeTopLevelItem(i)
                return l
            if l.childCount()>0:
                u = self._removeObj(l, obj)
                if u:
                    return u
        return None

    def _removeObj(self, l, obj):
        for i in range(l.childCount()):
            m = l.child(i)
            if m.obj == obj:
                l.takeChild(i)
                return m
            if m.childCount()>0:
                u = self._removeObj(m, obj)
                if u:
                    return u
        return None
    
    def toggleVis(self, item):
        print(item, item.prevstate)
        if item.prevstate==QtCore.Qt.Unchecked:
            item.setCheckState(0, QtCore.Qt.PartiallyChecked)

        item.prevstate = item.checkState(0)

        c = item.checkState(0)>0
        d = item.checkState(0)==1
        print('cs', item.checkState(0))

        update = False
        if c != item.obj.visible:
            item.obj.visible = c
            update = True
#            self.GLwidget.updateGL()
        if d != item.obj.transparent:
            item.obj.transparent = d
            update = True
#        if update:
#            self.GLwidget.updateGL()

    def updateParent(self, obj):
        bs = self.list_widget.blockSignals(True)
        parent = None #obj.parent
        s = self.list_widget.selectedItems()
        if parent == None:
            l = self.removeObject(obj)
            self.list_widget.addTopLevelItem(l)
            if l in s:
                self.setSelectedObj(obj)
        else:
            p_obj = self.getObjItem(parent)

            if not p_obj:
                print('Parent not found')
                return

            l = self.removeObject(obj)
            p_obj.addChild(l)
            if l in s:
                self.setSelectedObj(obj)
        self.list_widget.blockSignals(bs)

            

    def getObjItem(self, obj):
        for i in range(self.list_widget.topLevelItemCount()):
            l = self.list_widget.topLevelItem(i)
            if l.obj == obj:
                return l
            u = self._getObjItem(l, obj)
            if u:
                return u
        return None

    def _getObjItem(self, l, obj):
        for i in range(l.childCount()):
            m = l.child(i)
            if m.obj == obj:
                return m
            u = self._getObjItem(m, obj)
            if u:
                return u
        return None

    def setSelectedObj(self, obj):
        bs = self.list_widget.blockSignals(True)
        for i in range(self.list_widget.topLevelItemCount()):
            l = self.list_widget.topLevelItem(i)
            self._setSelectedObj(l, obj)
        self.list_widget.blockSignals(bs)
        self.update_label(obj)

    def _setSelectedObj(self, l, obj):
        if l.obj==obj:
            l.setSelected(True)
        else:
            l.setSelected(False)
        for i in range(l.childCount()):
            m = l.child(i)
            self._setSelectedObj(m, obj)

    def update_label(self, obj):
        self.label.setText(obj.info)
        self.label.setWordWrap(True)

    def getObj(self):
        return self.list_widget.currentItem().obj_name

    def itemSelectionChanged(self):
        l = self.list_widget.selectedItems()
        if l:
            pass
#            self.GLwidget.set_active_obj(l[0].obj_name)
        else:
            pass
#            self.GLwidget.set_active_obj(None)

def point_to_np(p):
    return np.array((p.x(), p.y()))

class StackItem(QtWidgets.QGraphicsPixmapItem):
    def __init__(self, stack):
        self.stack = stack
        self.pm = QtGui.QPixmap.fromImage(to_qimg(stack.maxproj))
        QtWidgets.QGraphicsPixmapItem.__init__(self, self.pm)

    def mousePressEvent(self, ev):
        print(ev)
        print(ev.scenePos())

    def mouseMoveEvent(self, ev):
        print(ev.scenePos(),ev.lastScenePos())
        delta = np.asarray(ev.scenePos()) - np.asarray(ev.lastScenePos())
        mode = self.rw.click_mode 
        if mode=='translate':
            self.stack.user_transform.matrix = np.array(((1,0,0,0),(0,1,0,delta.x()),(0,0,1,delta.y()),(0,0,0,1))).dot(self.stack.user_transform.matrix)
        else:
            c = self.boundingRect().center()
            c = np.array([0.0, c.x(), c.y(), 1.0])
            c = self.stack.user_transform.matrix.dot(self.stack.default_transform.matrix).dot(c)[1:3]
            d0 = point_to_np(ev.lastScenePos()) - c
            d1 = point_to_np(ev.scenePos()) - c
            if mode=='scale':
                s = la.norm(d1)/la.norm(d0)

                m = np.array(((1,0,0,0), (0,1,0,-c[0]),(0,0,1,-c[1]),(0,0,0,1)))
                r = np.array(((1,0,0,0), (0,1,0,c[0]),(0,0,1,c[1]),(0,0,0,1)))
                m = r.dot(np.array(((1,0,0,0),(0,s,0,0),(0,0,s,0),(0,0,0,1)))).dot(m)
                self.stack.user_transform.matrix = m.dot(self.stack.user_transform.matrix)
            else:
                theta = -math.asin((d1[1]*d0[0]-d1[0]*d0[1])/la.norm(d1)/la.norm(d0))
                m = np.array(((1,0,0,0), (0,1,0,-c[0]),(0,0,1,-c[1]),(0,0,0,1)))
                r = np.array(((1,0,0,0), (0,1,0,c[0]),(0,0,1,c[1]),(0,0,0,1)))

                s = math.sin(theta)
                c = math.cos(theta)
                rot = r.dot(np.array(((1,0,0,0),(0,c,s,0),(0,-s,c,0),(0,0,0,1.0)))).dot(m)
                self.stack.user_transform.matrix = rot.dot(self.stack.user_transform.matrix)
                    
                
        self.setTransform(QtGui.QTransform(*self.stack.get_transform2d()))
        
class RenderWidget(QtWidgets.QGraphicsView):
    def __init__(self):
        QtWidgets.QGraphicsView.__init__(self)#, self.scene)
        self.gscene = QtWidgets.QGraphicsScene(self)
        self.setScene(self.gscene)
        self.gscene.setBackgroundBrush(QtGui.QBrush(QtCore.Qt.black, QtCore.Qt.SolidPattern))
        self.click_mode = 'translate'
        self.stack_obj = {}


    def change_mode(self):
        print(self.sender())
        self.click_mode = self.sender().property('action')
        print(self.click_mode)
        
    def add_stack(self, stack):
        so = StackItem(stack)

#        QtWidgets.QGraphicsPixmapItem(QtGui.QPixmap.fromImage(to_qimg(stack.maxproj)))
        so.setTransform(QtGui.QTransform(*stack.get_transform2d()))
        so.setOpacity(0.7)
        so.rw = self
        self.stack_obj[stack] = so
        self.gscene.addItem(so)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)

        self.wc = WorldController()

        self.renderwidget = RenderWidget()
        self.setCentralWidget(self.renderwidget)
        self.tbar = QtWidgets.QToolBar()
        translate_button = QtWidgets.QPushButton('translate')
        translate_button.setProperty('action', 'translate')
        rotate_button = QtWidgets.QPushButton('rotate')
        rotate_button.setProperty('action', 'rotate')
        scale_button =  QtWidgets.QPushButton('scale')
        scale_button.setProperty('action', 'scale')
        self.tbar.addWidget(translate_button)
        self.tbar.addWidget(rotate_button)
        self.tbar.addWidget(scale_button)
        translate_button.clicked.connect(self.renderwidget.change_mode)
        rotate_button.clicked.connect(self.renderwidget.change_mode)
        scale_button.clicked.connect(self.renderwidget.change_mode)
        
        
        self.addToolBar(self.tbar)

        self.sbar = StackSideBar()

        self.wc.sbar = self.sbar
        self.wc.rw = self.renderwidget
        
        self.addDockWidget(QtCore.Qt.DockWidgetArea(2), self.sbar)
        self.make_menu()


    def action_load_stack(self):
        r = QtWidgets.QFileDialog.getOpenFileName(self, 'Stack', '.', 'TIFF files (*.tif)')
        if type(r)==tuple:
            r = r[0]
        if r:
            filename = str(r)
            self.wc.load_stack(filename)

    def quit_action(self):
        self.close()
        
    def make_menu(self):
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&File')
        load_action = QtWidgets.QAction('&Open Stack ...', self)
        load_action.setShortcut('Ctrl+O')
        quit_action = QtWidgets.QAction('&Quit', self)
        quit_action.setShortcut('Ctrl+Q')
        fileMenu.addAction(load_action)
        fileMenu.addAction(quit_action)
        load_action.triggered.connect(self.action_load_stack)
        #save_action.triggered.connect(self.action_write_label_stack)
        quit_action.triggered.connect(self.quit_action)
                                                                                
        

def main():

    parser = argparse.ArgumentParser(description='Generate surface mesh from stack')

    args = parser.parse_args()
    print(args)

    app = QtWidgets.QApplication(['Align-o-troN'])

    window = MainWindow()

    window.show()
    app.exec_()


if __name__=="__main__":
    main()
