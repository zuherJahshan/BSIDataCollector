import pandas as pd
from typing import List, Dict, Tuple
import os
from datetime import date
import gzip

class Sample(object):
    Property = str
    Value = str
    PropertyValueMap = Dict[Property, Value]
    def __init__(self,
                 property_value_map: PropertyValueMap) -> None:
        self.property_value_map = property_value_map
    
    def getProperties(self) -> List[Property]:
        return list(self.property_value_map.keys())
    
    def getValueOfProperty(self,
                           property: Property) -> Value:
        if not property in self.property_value_map:
            return None
        else:
            return self.property_value_map[property]
    
    # This function should change and become much more generic
    def downloadSample(self,
                        path_to_dir: str = None) -> None:
        if not 'fastq_ftp' in self.property_value_map or not 'sample_accession' in self.property_value_map:
            print('This Sample is not downloadable')
            return
        # Configure location and filename
        if path_to_dir == None:
            path_to_dir = './'
                    
        # download file
        os.makedirs(path_to_dir,
                    exist_ok = True)
        for filename, uri in self._getMissingFilesAndURIs(path_to_dir):
            filepath = path_to_dir + '/' + filename
            os.system('wget {} -O {}'.format(uri, filepath))
            os.system('gzip -d {}'.format(filepath))
    
    def getExistingLocalFiles(self,
                              path_to_dir):
        for filename, _ in self._getFilenamesAndURIs():
            unzipped_filename = '.'.join(filename.split('.')[:-1])
            filepath = path_to_dir + '/' + unzipped_filename
            if os.path.exists(filepath):
                yield filepath
    
    ######## Private ##########
    def _getMissingFilesAndURIs(self,
                                path_to_dir):
        # A sample will be available if all of its links are available
        for filename, uri in self._getFilenamesAndURIs():
            unzipped_filename = '.'.join(filename.split('.')[:-1])
            if not os.path.exists(path_to_dir + '/' + unzipped_filename):
                yield filename, uri
        
    
    def _getFilenamesAndURIs(self) -> str:
        URIs = self.property_value_map['fastq_ftp'].split(';')
        for uri in URIs:
            yield uri.split('/')[-1], uri
        
        
    

Patient = str
####################################
############ BSISample class #######
####################################
class BSISample(Sample):
    class SampleType(object):
        Blood = 'Blood'
        Capsule = 'Capsule'
        FMT = 'FMT'
    
    def __init__(self,
                 property_value_map: Sample.PropertyValueMap) -> None:
        super().__init__(property_value_map)
        self.id = self.property_value_map['sample_accession']
        self.type: BSISample.SampleType = None
        self.patient: BSISample.Patient = None
        self.date: date = 'NA'
        
    def getType(self) -> SampleType:
        if not self.type == None:
            pass
        else:
            if 'Blood' in self.property_value_map['experiment_alias']:
                self.type = BSISample.SampleType.Blood
            elif 'FMT' in self.property_value_map['experiment_alias']:
                self.type = BSISample.SampleType.FMT
            elif 'Capsule' in self.property_value_map['experiment_alias']:
                self.type = BSISample.SampleType.Capsule
        return self.type
    
    def getPatient(self) -> Patient:
        if not self.patient == None:
            return self.patient
        elif self.getType() == BSISample.SampleType.Capsule:
            self.patient = 'NA'
        else:
            self.patient = self.property_value_map['experiment_alias'].split('.')[0]
        return self.patient
    
    def getDate(self) -> date:
        if self.getType() == BSISample.SampleType.Capsule:
            self.date = 'NA'
        else:
            date_str = self.property_value_map['experiment_alias'].split('.')[1]
            self.date = date(year=int(date_str[0:4]),
                             month=int(date_str[4:6]),
                             day=int(date_str[6:]))
        return self.date
    
    def getID(self) -> str:
        return self.id
    
    def __str__(self):
        self.getType()
        self.getPatient()
        self.getDate()
        return "[\n\tSample ID:\t{}\n\tPatient:\t{}\n\tType:\t\t{}\n\tDate:\t\t{}\n]".format(self.id,
                                                                                             self.patient,
                                                                                             self.type,
                                                                                             self.date)
        
    def __repr__(self):
        return self.__str__()

Read = str
####################################
############ BSIDataCollector ######
####################################
class BSIDataCollector(object):
    BloodSample = BSISample
    FMTSample = BSISample
    CapsuleSample = BSISample
    PatientSamplesMap = Dict[Patient, Tuple[List[BloodSample], List[FMTSample]]]
    
    def __init__(self,
                 report_file: str = 'fastq_file_report.txt') -> None:
        
        # Load the file into a data frame
        df = pd.read_csv(report_file, sep='\t')
        
        # Build patient samples map, and capsules list
        self.capsule_samples: List[BSIDataCollector.CapsuleSample] = []
        self.patient_samples_map: BSIDataCollector.PatientSamplesMap = {}
        for row in df.iterrows():
            # build sample
            sample = BSISample(dict(row[1]))
            if sample.getType() == BSISample.SampleType.Capsule:
                self.capsule_samples.append(sample)
            else:
                self._addToPatientSampleMap(sample)
        
    def getPatients(self) -> List[Patient]:
        return list(self.patient_samples_map.keys())
    
    def getCapsules(self) -> List[CapsuleSample]:
        return self.capsule_samples
    
    def getBloodSamples(self,
                        patient) -> List[BSISample]:
        if not patient in self.patient_samples_map:
            return []
        else:
            return self.patient_samples_map[patient][0]
    
    def getFMTSamples(self,
                      patient) -> List[BSISample]:
        if not patient in self.patient_samples_map:
            return []
        else:
            return self.patient_samples_map[patient][1]
    
    def downloadSamples(self,
                        samples: List[BSISample]) -> None:
        for sample in samples:
            sample.downloadSample(path_to_dir=self._getSamplePath(sample))
            
    def getReads(self,
                 sample: BSISample) -> Read:
        sample.downloadSample(path_to_dir=self._getSamplePath(sample))
        for filepath in sample.getExistingLocalFiles(self._getSamplePath(sample)):
            with open(filepath) as f:
                read = [None, None, None]
                read_internal_idx = 0
                for line in f:
                    if read_internal_idx != 2:
                        read[read_internal_idx-int(read_internal_idx/3)] = line[:-1]
                    read_internal_idx += 1
                    read_internal_idx %= 4
                    if read_internal_idx == 0:
                        yield read
                
                    
                  
    ####################################
    ############ Private ###############
    ####################################
    def _addToPatientSampleMap(self,
                               sample: BSISample) -> None:
        assert sample.getType() != BSISample.SampleType.Capsule
        
        # If patient not found, add a new entry.
        if not sample.getPatient() in self.patient_samples_map:
            self.patient_samples_map.update({sample.getPatient(): ([],[])})
        
        # update the sample to the patient
        self.patient_samples_map[sample.getPatient()][int(sample.getType() == BSISample.SampleType.FMT)].append(sample)
        
    def _getSamplePath(self,
                       sample: BSISample) -> str:
        if sample.getType() == BSISample.SampleType.Capsule:
            return 'Data/Capsules/'
        elif sample.getType() == BSISample.SampleType.Blood:
            return 'Data/{}/Blood/'.format(sample.getPatient())
        elif sample.getType() == BSISample.SampleType.FMT:
            return 'Data/{}/FMT/'.format(sample.getPatient())