class index_labels:
    """ Class doc 
    
    A class structure for storing and retrieving row and index labels.
    
    """
    
    row_labels=[]
    column_labels=[]
    
    def __init__ (self):
        """ Class initialiser """
        pass
        
    def get_rows(self):
        """ Function doc """
        return self.row_labels

    def get_number_of_rows(self):
        """ Function doc """
        return len(self.row_labels)


    def get_columns(self):
        """ Function doc """
        return self.column_labels

    def get_number_of_columns(self):
        """ Function doc """
        return len(self.column_labels)
        

    def set_rows( self, new_row_labels ):
        """ Function doc """
        self.row_labels = new_row_labels

    def set_columns( self, new_row_labels ):
        """ Function doc """
        self.column_labels = new_row_labels
