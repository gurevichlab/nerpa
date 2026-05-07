pub struct SplitWhen<I, F>
where
    I: Iterator,
{
    iter: I,
    pred: F,
    pending: Option<I::Item>,
}

pub fn split_when<I, F>(iterable: I, pred: F) -> SplitWhen<I::IntoIter, F>
where
    I: IntoIterator,
    F: FnMut(&I::Item, &I::Item) -> bool,
{
    SplitWhen {
        iter: iterable.into_iter(),
        pred,
        pending: None,
    }
}

impl<I, F> Iterator for SplitWhen<I, F>
where
    I: Iterator,
    F: FnMut(&I::Item, &I::Item) -> bool,
{
    type Item = Vec<I::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        let first = match self.pending.take().or_else(|| self.iter.next()) {
            Some(x) => x,
            None => return None,
        };

        let mut chunk = vec![first];

        while let Some(item) = self.iter.next() {
            let prev = chunk.last().unwrap();
            if (self.pred)(prev, &item) {
                self.pending = Some(item); // start next chunk with this item
                break;
            }
            chunk.push(item);
        }

        Some(chunk)
    }
}
