################################################################################
#################     ImpulseDE2 object type unions     ########################
################################################################################

# Define class unions for slots
setClassUnion('numericORNULL', members = c('numeric', 'NULL'))
setClassUnion('listORNULL', members = c('list', 'NULL'))
setClassUnion('data.frameORNULL', members = c('data.frame', 'NULL'))