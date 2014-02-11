mkdir obsolete
rm *~
mv native* obsolete/.
mv junk* obsolete/.
mv modify* obsolete/.
mv tester* obsolete/.
