# -*- coding: utf-8 -*-

"""
***************************************************************************
    RandomPointsAlongEachLine.py
    ---------------------
    Date                 : February 2020
    Copyright            : (C) 2020 by Håvard Tveite
    Email                : havard dot tveite at nmbu dot no
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Håvard Tveite'
__date__ = 'February 2020'
__copyright__ = '(C) 2020, Håvard Tveite'

import random

from qgis.PyQt.QtCore import QVariant
from qgis.core import (QgsField,
                       QgsFeatureSink,
                       QgsFeature,
                       QgsFields,
                       QgsGeometry,
                       QgsPointXY,
                       QgsWkbTypes,
                       QgsSpatialIndex,
                       QgsFeatureRequest,
                       QgsDistanceArea,
                       QgsProject,
                       QgsProcessing,
                       QgsProcessingException,
                       QgsProcessingParameterDistance,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingOutputNumber,
                       QgsProcessingParameterDefinition)

from processing.algs.qgis.QgisAlgorithm import QgisAlgorithm
from processing.tools import vector


class HTRandomPointsAlongEachLine(QgisAlgorithm):

    INPUT = 'INPUT'
    POINTS_NUMBER = 'POINTS_NUMBER'
    MIN_DISTANCE = 'MIN_DISTANCE'
    OUTPUT = 'OUTPUT'
    MAXTRIESPERPOINT = 'MAX_TRIES_PER_POINT'
    OUTPUT_POINTS = 'OUTPUT_POINTS'

    def group(self):
        return self.tr('Vector creation')

    def groupId(self):
        return 'vectorcreation'

    def __init__(self):
        super().__init__()

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT,
                                                              self.tr('Input layer'),
                                                              [QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterNumber(self.POINTS_NUMBER,
                                                       self.tr('Number of points on each line'),
                                                       QgsProcessingParameterNumber.Integer,
                                                       1, False, 1, 1000000000))
        self.addParameter(QgsProcessingParameterDistance(self.MIN_DISTANCE,
                                                         self.tr('Minimum distance between points'),
                                                         0, self.INPUT, False, 0, 1000000000))
        self.addParameter(QgsProcessingParameterNumber(self.MAXTRIESPERPOINT,
                                                       self.tr('Maximum number of attempts per point'),
                                                       QgsProcessingParameterNumber.Integer,
                                                       10, False, 1, 1000))
        self.addParameter(QgsProcessingParameterFeatureSink(self.OUTPUT,
                                                            self.tr('Random points'),
                                                            type=QgsProcessing.TypeVectorPoint))

        self.addOutput(QgsProcessingOutputNumber(
                         self.OUTPUT_POINTS,
                         self.tr('Number of point generated')
                       ))



    def name(self):
        return 'htrandompointsalongeachline'

    def displayName(self):
        return self.tr('HT Random points along each line')

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT))

        pointCount = self.parameterAsDouble(parameters, self.POINTS_NUMBER, context)
        minDistance = self.parameterAsDouble(parameters, self.MIN_DISTANCE, context)
        MaxTriesPerPoint = self.parameterAsDouble(parameters, self.MAXTRIESPERPOINT, context)
        fields = QgsFields()
        fields.append(QgsField('id', QVariant.Int, '', 10, 0))

        (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT, context,
                                               fields, QgsWkbTypes.Point, source.sourceCrs(), QgsFeatureSink.RegeneratePrimaryKey)
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT))
      
        totNPoints = 0  # The total number of points generated
        featureCount = source.featureCount()
        total = 100.0 / (pointCount * featureCount) if pointCount else 1
        random.seed()

        index = QgsSpatialIndex()
        points = dict()

        da = QgsDistanceArea()
        da.setSourceCrs(source.sourceCrs(), context.transformContext())
        da.setEllipsoid(context.project().ellipsoid())

        maxIterations = pointCount * MaxTriesPerPoint
        for f in source.getFeatures():
            lineGeoms = []
            lineCount = 0
            fGeom = f.geometry()
            feedback.pushInfo('fGeom: ' + str(fGeom))
            totLineLength = da.measureLength(fGeom)
            feedback.pushInfo('fGeom totLineLength: ' + str(totLineLength))
            # Explode multi part
            if fGeom.isMultipart():
                for aLine in fGeom.asMultiPolyline():
                    lineGeoms.append(aLine)
                #lines = fGeom.asMultiPolyline()
                # pick random line
                #lineId = random.randint(0, len(lines) - 1)
                #vertices = lines[lineId]
            else:
                lineGeoms.append(fGeom.asPolyline())
                #vertices = fGeom.asPolyline()
            feedback.pushInfo('lineGeoms: ' + str(lineGeoms))


            # Generate points on the line geometry / geometries
            nPoints = 0
            nIterations = 0
            while nIterations < maxIterations and nPoints < pointCount:
                if feedback.isCanceled():
                    break
                #feedback.pushInfo('nIterations: ' + str(nIterations))
                # Get the random "position" for this point
                randomLength = random.random() * totLineLength
                feedback.pushInfo('randomLength: ' + str(randomLength))
                currLength = 0
                prefLength = 0
                # Go through the parts
                for l in lineGeoms:
                    if feedback.isCanceled():
                        break
                    currGeom = QgsGeometry.fromPolylineXY(l)
                    #lineLength = da.measureLength(QgsGeometry.fromPolylineXY(l))
                    lineLength = da.measureLength(currGeom)
                    prevLength = currLength
                    currLength += lineLength
                    feedback.pushInfo('l lineLength: ' + str(lineLength) + ' currLength: ' + str(currLength))
                    vertices = l
                    # Skip if this is not the "selected" part
                    if currLength < randomLength:
                        continue
                    #randomLength -= currLength
                    #vertices = QgsGeometry.fromPolylineXY(l)
                    feedback.pushInfo('l/vertices: ' + str(vertices))

                    #randomLength = random.random() * lineLength
                    #distanceToVertex(vid)

                    # find the segment for the new point
                    # and calculate the offset (remainDistance) on that segment
                    remainDist = randomLength - prevLength
                    feedback.pushInfo('remainDist1: ' + str(remainDist))
                    if len(vertices) == 2:
                        vid = 0
                        #remainDist = randomLength - currLength
                    else:
                        vid = 0
                        #while (fGeom.distanceToVertex(vid)) < randomLength:
                        #while (currGeom.distanceToVertex(vid)) < randomLength:
                        currDist = currGeom.distanceToVertex(vid)
                        prevDist = currDist
                        while currDist < remainDist and vid < len(vertices):
                            vid += 1
                            prevDist = currDist
                            currDist = currGeom.distanceToVertex(vid)
                            feedback.pushInfo('currdist: ' + str(currDist) + ' vid: ' + str(vid))
                        if vid == len(vertices):
                            feedback.pushInfo('**** vid = len(vertices)! ****')
                        vid -= 1
                        feedback.pushInfo('currdist2: ' + str(currDist) + ' vid: ' + str(vid))
                        remainDist = remainDist - prevDist
                    feedback.pushInfo('remainDist2: ' + str(remainDist))
                    if remainDist <= 0:
                        continue
                    startPoint = vertices[vid]
                    endPoint = vertices[vid + 1]
                    length = da.measureLine(startPoint, endPoint)
                    # if remainDist > minDistance:
                    d = remainDist / (length - remainDist)
                    rx = (startPoint.x() + d * endPoint.x()) / (1 + d)
                    ry = (startPoint.y() + d * endPoint.y()) / (1 + d)

                    # generate random point
                    p = QgsPointXY(rx, ry)
                    geom = QgsGeometry.fromPointXY(p)
                    if vector.checkMinDistance(p, index, minDistance, points):
                        f = QgsFeature(totNPoints)
                        f.initAttributes(1)
                        f.setFields(fields)
                        f.setAttribute('id', totNPoints)
                        f.setGeometry(geom)
                        sink.addFeature(f, QgsFeatureSink.FastInsert)
                        index.addFeature(f)
                        points[nPoints] = p
                        nPoints += 1
                        totNPoints += 1
                        feedback.setProgress(int(totNPoints * total))
                    break
                nIterations += 1

            if nPoints < pointCount:
                #feedback.pushInfo(self.tr('Could not generate requested number of random points. '
                feedback.reportError(self.tr('Could not generate requested number of random points. '
                                             'Maximum number of attempts exceeded.'), False)

        return {self.OUTPUT: dest_id, self.OUTPUT_POINTS: nPoints}
